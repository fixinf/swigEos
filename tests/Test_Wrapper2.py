import unittest
from Wrapper2 import Model
import Models2
import numpy as np
from os.path import join
import matplotlib.pyplot as plt
from numpy.testing import assert_allclose, assert_almost_equal

class TestClassicMods(unittest.TestCase):
    def setUp(self):
        self.ref_dir = '/home/const/Dropbox/Documents/For DN/Very final/data'
        self.KV = Models2.KVOR()
        self.W = Models2.waleckaMatsui()
#         self.W_params = np.loadtxt(join(self.ref_dir, 'waleckaMatsui','props.dat'),
#                                    skiprows=1, usecols=[1],delimiter='   ')
#         self.KVOR_params = np.loadtxt(join(self.ref_dir, 'KVOR','props.dat'),
#                                    skiprows=1, usecols=[1])
        self.W_eos = np.loadtxt(join(self.ref_dir, 'waleckaMatsui','eos.dat'),
                                skiprows=1)
        self.KV_eos = np.loadtxt(join(self.ref_dir, 'KVOR','eos.dat'),
                                skiprows=1)
#         self.W_eos = np.loadtxt(join(self.ref_dir, 'waleckaMatsui','eos.dat'),
#                                 skiprows=1)
        self.KV_hyper = np.loadtxt(join(self.ref_dir, 'KVOR','hyper.dat'),
                                skiprows=1)
        self.KV_mcrust = np.loadtxt(join(self.ref_dir, 'KVOR','masses_crust.dat'),
                                skiprows=1)
        self.KV_mcrust_hyper = np.loadtxt(join(self.ref_dir, 'KVOR','mass_hyper.dat'),
                                skiprows=1)
    def eos_compare(self, n, m, t):
        np.testing.assert_allclose(m.sym.Ebind(n), t[:,1], rtol=1e-5, atol=0)
#         plt.plot(n, m.sym.P(n), n, t[:,2])
#         plt.show()
#         print m.sym.P(n) - t[:,2]
        np.testing.assert_allclose(t[:,2], m.sym.P(n),  rtol=1e-5, atol=1e-4, 
                                   verbose=1)
        np.testing.assert_allclose(m.neutr.Ebind(n), t[:,8], rtol=1e-5, atol=0)
        np.testing.assert_allclose(t[:,9], m.neutr.P(n),  rtol=1e-5, atol=1e-4,
                                   verbose=1)


    def testWaleckaEos(self): 
        pass#self.eos_compare(self.W_eos[:, 0] * self.W.n0, self.W, self.W_eos)      
        
    def testWaleckaParams(self):
        pass

    def testKVOR(self):
        return
        self.eos_compare(self.KV_eos[:, 0] * self.KV.n0, self.KV, self.KV_eos)
    
    def testUopt(self):
        pass
    
    def nucl_compare(self, m, t):
        nmin = t[0, 0] * m.n0
        nmax = t[-1,0] * m.n0
        npoints = t.shape[0] - 1
        # m.nucl.reset(nmin=nmin, nmax=nmax, npoints=npoints)
        # m.nucl.reset(nrange=t[:,0]*m.n0)

        E, P, n = m.nucl.EPN(nrange=t[:,0]*m.n0)

        print(m.nucl.nrange)
        print(t[:, 0] * m.n0)
        np.testing.assert_allclose(E, t[:,4], rtol=1e-5)
        np.testing.assert_allclose(P, t[:,5], rtol=1e-5)
        np.testing.assert_allclose(m.nucl.concentrations()[:, 1], t[:,6],
                                   rtol=1e-5, atol=1e-4)
        np.testing.assert_allclose(m.nucl.rho[:, 0], t[:,7],
                                   rtol=1e-5, atol=1e-4)

    
    def hyper_compare(self, m, t):
        nmin = t[0, 0] * m.n0
        nmax = t[-1,0] * m.n0
        npoints = t.shape[0] - 1
        decimal = 4
        # m.hyper.reset(nmin=nmin, nmax=nmax, npoints=npoints)

        E, P, n = m.hyper.EPN(nrange=t[:,0]*m.n0)
        print(m.hyper.nrange)
        print(t[:, 0] * m.n0)
        assert_almost_equal(m.hyper.nrange, t[:, 0]*m.n0)
        assert_almost_equal(E, t[:,1], decimal=decimal)
        assert_almost_equal(P, t[:,2], decimal=decimal)
        assert_almost_equal(m.hyper.concentrations(), t[:,3:11],
                                   decimal=decimal)
        assert_almost_equal(m.hyper.rho[:, 0], t[:,11],
                                   decimal=decimal)
        
    def mcrust_compare(self, m=Models2.KVOR(), t=None):
        n_stars = t[:, 0] * m.n0
        kwargs = {'rtol':1e-4, 'atol':1e-3}
        out = np.array(m.nucl.stars_crust(n_stars=n_stars, inter='cubic')).transpose()
        assert_allclose(out[:, 0], t[:, 0]*m.n0, **kwargs)
        assert_almost_equal(out[:, 1], t[:, 1], decimal=3)
        assert_almost_equal(out[:, 2], t[:, 2], decimal=3)
        assert_allclose(out[:, 4], t[:, 4], **kwargs)

    def mcrust_compare_hyper(self, m=Models2.KVOR(), t=None):
        # return
        n_stars = t[:, 0] * m.n0
        kwargs = {'rtol':1e-4, 'atol':1e-3}
        out = np.array(m.hyper.stars_crust(n_stars=n_stars, inter='cubic')).transpose()
        assert_allclose(out[:, 0], t[:, 0]*m.n0, **kwargs)
        assert_almost_equal(out[:, 1], t[:, 1], decimal=3)
        assert_almost_equal(out[:, 2], t[:, 2], decimal=3)
        assert_allclose(out[:, 4], t[:, 4], **kwargs)

    def testWnucl(self):
        return
        self.nucl_compare(self.W, self.W_eos)

    def testKVORnucl(self):
        # return
        self.nucl_compare(self.KV, self.KV_eos)

    def testKVORhyper(self):
        return
        self.hyper_compare(self.KV, self.KV_hyper)
        self.KV.hyper.dumpMeff()

    def testKVORmcrust(self):
        return
        self.mcrust_compare(self.KV, self.KV_mcrust)

    def testKVORmhyper(self):
        return
        self.mcrust_compare_hyper(self.KV, self.KV_mcrust_hyper)

    def testKVORdumpEOS(self):
        return
        self.KV.dumpEos()
        self.KV.nucl.dumpMeff()

    def testKVORdumpHyperEOS(self):
        return
        self.KV.hyper_phi_sigma.dumpEos()
        self.KV.hyper_phi_sigma.dumpMeff()




if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
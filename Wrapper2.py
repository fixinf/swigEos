import eosWrap as eos
import numpy as np
from numpy import array as arr, pi, sqrt
import matplotlib.pyplot as plt
import os
from os.path import join
from tabulate import tabulate
from scipy import interpolate
from scipy.misc import derivative
from multiprocessing import Queue, Process
from Wrapper import Wrapper as Wrapper_old
import six

BASEFOLDER = '/home/const/wr2_test'

class Wrapper(object):
    def __init__(self, Ctype):
        self.Ctype = Ctype
        self.C = Ctype()
        self.foldername = os.path.join(BASEFOLDER, Ctype.__name__.strip('_'))
        if not os.path.exists(self.foldername):
            os.makedirs(self.foldername)
        self.m_pi = 135.
        self.m_e = 0.5 / self.m_pi
        self.m_mu = 105.0 / self.m_pi
        self.n0 = self.C.n0
        self.set = False
        self.driver_set = False
        self.folder = None
        self.n_baryon = 2
        self.verbose = False
        self.mpi4_2_mevfm3 = self.m_pi * (self.m_pi / 197.33) ** 3
        self.nmin = 0.
        self.nmax = 10*self.n0
        self.npoints = 1000
        self.nrange = np.linspace(self.nmin, self.nmax, self.npoints,
                                  endpoint=False)
        #example of setting filenames for a descendant of Wrapper
        self.filenames = {'mass_nocrust': None, 'eos': None,
                          'mass_crust' : None, 'meff':None,
                          'eta' : None, 'vs' : None,
                          'parts':None}

        self.part_names = ['n', 'p', 'Lambda', 'Sigma-',
                           'Sigma0', 'Sigma+', 'Xi-', 'Xi0', 'e', 'mu']


    def stepE(self, n, last, f, lenLast, iter, C, que):
        rho = eos.stepE(n, last, f, len(last), iter, C)
        que.put(rho)

    def reset(self, iterations=30, timeout=5):
        """Calculates the EoS stored in the wrapper. Has to be overrided in
        Symm, Neutron. 
        
        Calculates the particle composition, scalar field, energy density, 
        pressure, ... for the given parameter set self.C. Wrapper.set is 
        set to True after s succesful calculation. 
        
        Note: \sigma^* is not supported in the current version.

        """

        print "Resetting " + self.__repr__() + " for model " + self.Ctype.__name__.strip('_')
        init = arr([0. for i in range(self.n_baryon - 1)])

        rho = []

        mu_e = []
        f = arr([0.])  # sp omitted
        for i, _n in enumerate(self.nrange):

            if timeout is None:
                init = eos.stepE(_n, init, f, len(init), iterations, self.C)
            else:
                queue = Queue()
                p = Process(target=self.stepE, args=(_n, init, f, len(init),
                                                  iterations, self.C, queue))
                p.start()
                p.join(timeout)
                if p.is_alive():
                    p.terminate()
                    print "timeout reached"
                    self.rho = np.ascontiguousarray(rho[:])
                    self.nrange = self.nrange[: self.rho.shape[0]]
                    _E = map(lambda z: eos.E(z, self.C), self.rho)
                    self._E = np.ascontiguousarray(np.array(_E[:]))
                    self._P = np.ascontiguousarray(self.P_chem(self.rho))
                    self.set=1
                    return
                init = queue.get(timeout=None)
                # print init

            if i % (len(self.nrange) / 20) == 0:
                print '.',

            rho.append(init.copy())
            rho[i] = np.insert(rho[i], 0, _n - np.sum(init))
            f = eos.f_eq(rho[i], f, 1, self.C)  # sp omitted
            if self.verbose:
                pass  # TODO: some verbose output
            # rho contains all scalar fields as well
            rho[i] = np.insert(rho[i], 0, f)  #and again sp omitted
            mu_e.append(eos.mu(rho[i], 1, self.C) -
                        eos.mu(rho[i], 2, self.C))

        self.rho = np.ascontiguousarray(arr(rho))
        self.mu_e = np.ascontiguousarray(arr(mu_e))
        eparts = []
        _E = []
        for z in self.rho:
            epart = np.zeros((9), dtype='float')
            _E.append(eos.E(z, self.C, epart))
            eparts.append(epart)
        self._E = np.array(_E)
        self.Eparts = arr(eparts)
        # self._E = np.array(map(lambda z: eos.E(z, self.C), self.rho))
        self._E = np.ascontiguousarray(self._E)
        self._P = np.ascontiguousarray(self.P_chem(self.rho))
        self.set = True

    def EPN(self, nrange=None):
        if self.check(nrange=nrange):
            return self._E, self._P, self.nrange
        else:
            raise

    def check(self, nrange=None):  # TODO Rewrite this hell
        """Returns if the wrapper is set or not. If nrange is None and nmin, 
        nmax and npoints are all present, constructs nrange from nmin, nmax,
        npoints. Otherwise, nrange's priority is higher that other three 
        arguments' """
        if nrange is not None:
            if nrange.shape == self.nrange.shape:
                if np.allclose(nrange, self.nrange) and self.set:
                    return True
                else:
                    self.nrange = nrange
                    self.reset()
            else:
                self.nrange = nrange
                self.reset()
        else:
            if not self.set:
                self.reset()

        return self.set

    def solve(self, f0=0.195, E0=-16., K0=275., J0=32., mu_scale=1.):
        res = eos.solve(self.n0, E0, f0, K0, J0, self.C, iter, mu_scale)
        if self.verbose:
            pass  # TODO: some verbose output
        if mu_scale < 50000 and res == 5:
            self.solve(f0, E0, K0, J0, iter, 5 * mu_scale)

    def dumpMasses(self, nmin=0.4, nmax=4., npoints=100, write=False):
        if write:
            if not os.path.exists(self.foldername):
                os.makedirs(self.foldername)
        NPOINTS_EOS = 1000
        if self.check(nmin=0., nmax=nmax, npoints=NPOINTS_EOS):
            if not self.driverSet:
                self.setDriver()

        n, m, r, mb = self.stars(npoints=100, nmax=10. * self.n0)
        if write:
            table = np.array([n / self.n0, m, r, mb]).transpose()
            f = open(os.path.join(self.foldername, 'masses.dat'), 'w')
            tab = tabulate(table, ['n/n_0',
                                   'M [M_{sun}]',
                                   'R [km]',
                                   'M_B [M_sun]'],
                           tablefmt='plain')
            f.write(tab)
            f.close()
        return n, m, r


    def dumpMassesCrust(self, nmin=0.4, nmax=5., npoints=100, write=True, fname=None, ret_frac=True):
        self.check()
        out = self.stars_crust(nmin=nmin, nmax=nmax, npoints=npoints, ret_frac=ret_frac)
        out[0] /= self.n0
        if write:
            tab = arr(out).transpose()
            names = ['n/n0', 'M[M_sun]', 'R[km]', 'Mg_rect[M_sun]', 'Mg_trap[M_sun]']
            if ret_frac:
                names += self.part_names[1:self.n_baryon]
            table = tabulate(tab, names, tablefmt='plain')
            if fname is None:
                fname = self.filenames['mass_crust']
            with open(join(self.foldername, fname), 'w') as f:
                f.write(table)
        return out


    def stars_crust(self, ncut_crust=0.45, ncut_eos=0.6, inter='cubic',
                    n_stars=None, nmin=.6, nmax=5.0, npoints=50,
                    crust="crust.dat", show=False, crustName=None,
                    ret_frac=False, fasthyp=False, neutron=0):
        neos = self.npoints
        if n_stars is not None:
            nmax = n_stars[-1]  # TODO too dirty workaround
        self.check(nrange=self.nrange)
        E, P, N = self.EPN(self.nrange)
        E *= self.m_pi ** 4
        P *= self.m_pi ** 4
        N = N / self.n0
        e = []
        p = []
        n = []
        with open("/home/const/workspace/swigEosWrapper/crust.dat", 'r') as f:
            for line in f:
                # print line
                _e, _p, _n = line.split()
                if float(_n) < ncut_crust:
                    e.append(float(_e))
                    p.append(float(_p))
                    n.append(float(_n))
        n_eos = 5
        i_n_eos = np.argmin(abs(N - [ncut_eos for i in N]))

        plist = np.append(p[:], P[i_n_eos:(i_n_eos + n_eos)])
        elist = np.append(e[:], E[i_n_eos:(i_n_eos + n_eos)])
        nlist = np.append(n[:], N[i_n_eos:(i_n_eos + n_eos)])

        print i_n_eos
        print P[i_n_eos:(i_n_eos + n_eos)]
        print N[i_n_eos:(i_n_eos + n_eos)]
        # exit()
        iP = interpolate.interp1d(nlist, plist, kind=inter)
        iE = interpolate.interp1d(nlist, elist, kind=inter)

        gamma = 1. / 4.
        iN = np.linspace(1e-10 ** gamma, ncut_eos ** gamma, 10000)
        iN = iN ** (1. / gamma)

        crust_p = np.array(map(iP, iN))
        crust_e = np.array(map(iE, iN))

        #         finalE = np.append(crust_e, E[i_n_eos+n_eos:])
        #         finalP = np.append(crust_p, P[i_n_eos + n_eos:])
        #         finalN = np.append(iN, N[i_n_eos+n_eos:])

        finalE = np.append(crust_e, E[i_n_eos + n_eos:]) / self.m_pi ** 4
        finalP = np.append(crust_p, P[i_n_eos + n_eos:]) / self.m_pi ** 4
        finalN = np.append(iN, N[i_n_eos + n_eos:])

        if show:
            lines = plt.plot(finalN, finalP * self.m_pi ** 4, n, p, N, P)
            plt.xlabel(r'$n/n_0$', fontsize=18)
            plt.xlim([0, 2])
            plt.ylabel(r'$P \, [MeV^4] $', fontsize=18)
            plt.legend(lines, ['interpolated', 'crust', 'RMF'], loc=0)
            plt.show()

            lines = plt.plot(finalN, finalE * self.m_pi ** 4, n, e, N, E)
            plt.xlabel(r'$n/n_0$', fontsize=18)
            plt.xlim([0, 2])
            plt.ylabel(r'$E \, [MeV^4] $', fontsize=18)
            plt.legend(lines, ['interpolated', 'crust', 'RMF'], loc=0)
            plt.show()

        if crustName is not None:
            tab = np.array([finalN, finalE, finalP]).transpose()
            table = tabulate(tab, ['n/n0', 'E', 'P'], tablefmt='plain')
            with open(crustName, 'w') as f:
                f.write(table)

        self.dr = eos.KVDriver()
        self.dr.set(finalE, finalP, finalN * self.n0)
        if n_stars is None:
            n_stars = np.linspace(nmin, nmax, npoints)

        #interpolating particle fractions (except n)
        conc = self.concentrations()
        inter_hyp = [interpolate.interp1d(self.nrange, conc[:, i])
                         for i in range(1, self.n_baryon)]

        #         plt.plot(self.n, map(lambda z: [f(z) for f in inter_hyp], self.n))
        #         plt.show()
        MR = []
        Mgrav = []
        str_fracs = []
        for _n in n_stars:
            MR.append(eos.star_crust2(_n, 3, self.dr, 1e-11))
            lastN = self.dr.getLastN(self.dr.nSize)[:-1]
            lastR = self.dr.getLastR(self.dr.nSize)[:-1]
            lastM = self.dr.getLastM(self.dr.nSize)[:-1]
            dx = lastR[1] - lastR[0]

            grav_mult = []
            for i, r in enumerate(lastR):
                grav_mult.append(r ** 2 / sqrt(1 - 2 * 1.4677 * lastM[i] / r))

            grav_mult = np.array(grav_mult)
            res = np.multiply(lastN, grav_mult)
            integral = np.trapz(res, dx=dx)
            if ret_frac:
                hyper_N = []
                hyper_NR = []
                for f in inter_hyp:
                    hyper_N.append(lastN * map(f, lastN))
                str_frac = [np.trapz(np.multiply(hyper_N[i], grav_mult), dx=dx) / (3 * integral)
                            for i in range(0, self.n_baryon-1)]

                str_fracs.append(str_frac)

            Mgrav.append((0.0004898007281478712) * integral)

        MR = np.array(MR)
        Mgrav = np.array(Mgrav)
        str_fracs = np.array(str_fracs)

        if ret_frac:
            return [n_stars, MR[:, 0], MR[:, 1], 931.5 / self.m_pi * MR[:, 2],
                    931.5 / self.m_pi * Mgrav] + np.array(str_fracs).transpose().tolist()
        else:
            return (n_stars, MR[:, 0], MR[:, 1], 931.5 / self.m_pi * MR[:, 2],
                    931.5 / self.m_pi * Mgrav)


    def dumpScalings(self):
        E, f = self.E(self.nrange, ret_f=1)
        C = self.C
        tab = arr(map(lambda z: [z, C.eta_s(z), C.eta_o(z), C.eta_r(z), C.phi_n(0, z), C.U(z)], f))
        table = tabulate(np.insert(tab, 0, self.nrange/self.n0, axis=1), ['n/n0', 'f', 'eta_s', 'eta_o',
                                                                          'eta_r', 'phi_n', 'U'])
        with open(join(self.foldername, self.filenames['eta']), 'w') as f:
            f.write(table)

    """A stub. Descendants should provide their energy density through
    this method."""

    def P(self, n):
        pass

    """Returns pressure array for the given compositions
      n = {f, n_n, n_p, ...}. DO NOT OVERRIDE."""

    def P_chem(self, nlist, leptons=True):
        res = []
        sp = 1
        for n in nlist:
            n = arr(n)
            mu_e = eos.mu(n, sp, self.C) - eos.mu(n, sp + 1, self.C)
            n_e = 0.
            n_mu = 0.
            if (mu_e ** 2 - self.m_e ** 2 > 0):
                n_e += (mu_e ** 2 - self.m_e ** 2) ** (1.5) / (3 * pi ** 2)

            if (mu_e ** 2 - self.m_mu ** 2 > 0):
                n_mu += (mu_e ** 2 - self.m_mu ** 2) ** (1.5) / (3 * pi ** 2)

            if leptons:
                E = eos.E(n, self.C)
            else:
                E = eos._E(n, self.C)

            sum = 0.

            for i, r in enumerate(n[sp:]):
                sum += r * eos.mu(n, i + sp, self.C)

            if leptons:
                sum += n_e * mu_e + n_mu * mu_e

            res.append(sum - E)
        return np.array(res) * self.mpi4_2_mevfm3


    def E_gen(self, nlist, ret_f=False, f=0.):
        """Provides the energy density for given particle concentrations [n] 
        without leptons. DO NOT OVERRIDE.
        
        Input:
        n -- np.array; input baryon densities
        ret_f -- bool : do or not return f(n)
        f -- double : initial f(n) value
        Output:
        if ret_f = False --> np.array: E([n]) [m_\pi^4]
        else --> np.array, np.array: E([n]) [m_\pi^4],
                                     f([n])
        """
        flist = []
        elist = []
        eparts = []
        for n in nlist:
            epart = np.zeros(9, dtype='float')
            f = eos.f_eq(arr(n), arr([f]), 1, self.C)[0]
            flist.append(f)
            elist.append(eos._E(np.insert(n, 0, f), self.C, epart))
            eparts.append(epart)

        self.Eparts = arr(eparts)
        if ret_f:
            return [arr(elist), arr(flist)]
        else:
            return [arr(elist)]


    def E(self, nlist, ret_f=False, f=0.):
        """A stub. Descendants should provide their energy density through
        this method."""
        pass

    def Ebind(self, nrange, ret_f=False):
        out = self.E(nrange, ret_f=ret_f)
        res = np.nan_to_num(self.m_pi * (out[0] / nrange - self.C.M[0]))
        if ret_f:
            return res, out[1]
        else:
            return res

    def dumpLP(self):
        pass

    def concentrations(self):
        if not self.set:
            print 'Concentrations: wrapper is not set! Pass.'
            return

        rho = []
        sp = 1 + self.C.sprime
        for _r in self.rho:
            if np.sum(_r[sp:]) > 0:
                rho.append(_r[sp:] / np.sum(_r[sp:]))
            else:
                rho.append([1.] + [0. for i in _r[sp + 1:]])

        return np.array(rho)

    def dumpVs(self, nrange=None):
        E, P, n = self.EPN(nrange=nrange)
        vs = np.gradient(P) / np.gradient(E) / self.mpi4_2_mevfm3
        with open(join(self.foldername, self.filenames['vs']), 'w') as f:
            f.write(tabulate(arr([n/self.n0, vs]).transpose(), ['n/n0', 'v_s^2'], tablefmt='plain'))

    def dumpParts(self):
        self.check()
        if self.filenames['parts']:
            table = tabulate(np.insert(self.Eparts, 0, self.nrange/self.n0, axis=1),
                             ['n/n0', 'f', 'U(f)', 'K_n', 'K_p', 'om', 'rho',
                              'phi', 'e', 'mu'], tablefmt='plain')
            with open(join(self.foldername, self.filenames['parts']), 'w') as f:
                f.write(table)



class Nucleon(Wrapper):
    def __init__(self, C):
        super(Nucleon, self).__init__(C)
        self.C.Hyper = 0
        self.filenames['meff'] = 'meff.dat'
        self.filenames['mass_crust'] = 'masses_crust.dat'
        self.filenames['eta'] = 'eta_ns.dat'
        self.filenames['vs'] = 'vs_ns.dat'
        self.filenames['parts'] = 'parts_ns.dat'


    def E(self, nrange, ret_f=False, f=0.):
        E, P, n = self.EPN(nrange=nrange)
        if ret_f:
            flist = self.rho[:, 0]
            return E, flist
        else:
            return E

    def dumpEos(self):
        pass

    def lepton_concentrations(self, ret_du=False):
        """Returns lepton concentrations for self.nrange. """
        self.check()
        ne_list = []
        nmu_list = []
        du_list = []
        for mu in self.mu_e:
            n_e = 0
            n_mu = 0
            if mu > self.m_e:
                n_e += (mu**2 - self.m_e**2)**(1.5)/(3*pi**2)
            if mu > self.m_mu:
                n_mu += (mu**2 - self.m_mu**2)**(1.5)/(3*pi**2)

            ne_list.append(n_e)
            nmu_list.append(n_mu)

            if ret_du:
                if n_e + n_mu > 0:
                    du_list.append(1./(1 + (1 + (n_e/(n_e + n_mu))**(1./3))**3))
                else:
                    du_list.append(0.11)
        ne_list = arr(ne_list)
        nmu_list = arr(nmu_list)
        du_list = arr(du_list)
        if ret_du:
            return [ne_list, nmu_list, du_list]
        else:
            return [ne_list, nmu_list]

    def dumpMeff(self):
        self.check()
        tab = map(lambda z: [self.C.phi_n(i, self.C.Xs(i, z)*z)
                  for i in range(self.n_baryon)], self.rho[:, 0])
        data = np.insert(arr(tab), 0, self.nrange/self.n0, axis=1)
        table = tabulate(data, ['n/n0']+[self.part_names[i] for i in range(self.n_baryon)],
                         tablefmt='plain')
        with open(join(self.foldername, self.filenames['meff']), 'w') as f:
            f.write(table)



class Sym(Nucleon):
    def __init__(self, C):
        super(Sym, self).__init__(C)
        self.C.Hyper = 0
        self.filenames.update(J='j_tilde.dat')
        self.filenames['eta'] = 'eta_sym.dat'
        self.filenames['vs'] = 'vs_sym.dat'
        self.filenames['parts'] = 'parts_sym.dat'


    def reset(self):
        self._E, = self.E(self.nrange)
        self._P = self.P(self.nrange)
        self.set = 1

    def E(self, nrange, ret_f=False, f=0.):
        nlist = [[n / 2, n / 2] for n in nrange]
        return Wrapper.E_gen(self, nlist, ret_f, f)

    def P(self, nrange):
        nlist = []
        f = 0.
        for _n in nrange:
            f, = eos.f_eq(np.array([_n / 2, _n / 2]), np.array([f]), 1, self.C)
            nlist.append(np.array([f, _n / 2, _n / 2]))
        res = self.P_chem(nlist, leptons=False)
        return self.mpi4_2_mevfm3 * res

    def Jtilde(self, nrange=None):
        if nrange == None:
            nrange = self.nrange
        return arr(map(lambda z: eos.J(z, self.C), nrange))

    def dumpJ(self):
        tab = arr([self.nrange/self.n0, self.Jtilde(self.nrange)]).transpose()
        with open(join(self.foldername, self.filenames['J']), 'w') as f:
            f.write(tabulate(tab, ['n/n0', 'J\tilde[MeV]'], tablefmt='plain'))




class Neutr(Nucleon):
    def __init__(self, C):
        super(Neutr, self).__init__(C)
        self.C.Hyper = 0
        self.filenames['mass_crust'] = 'masses_crust_neutron.dat'
        self.filenames['mass_nocrust'] = 'masses_neutron.dat'
        self.filenames['meff'] = 'meff_ns.dat'
        self.filenames['eta'] = 'eta_n.dat'
        self.filenames['vs'] = 'vs_n.dat'
        self.filenames['parts'] = 'parts_neutr.dat'

    def reset(self):
        self._E, = self.E(self.nrange)
        self._P = self.P(self.nrange)
        self.set = 1

    def E(self, nrange, ret_f=False, f=0.):
        nlist = [[n, 0] for n in nrange]
        return Wrapper.E_gen(self, nlist, ret_f, f)

    def P(self, nrange):
        nlist = []
        f = 0.
        for _n in nrange:
            f, = eos.f_eq(np.array([_n, 0]), np.array([f]), 1, self.C)
            nlist.append(np.array([f, _n, 0]))
        res = self.P_chem(nlist, leptons=False)
        return self.mpi4_2_mevfm3 * res




class Hyperon(Nucleon):
    def __init__(self, C):
        super(Hyperon, self).__init__(C)
        self.C.Hyper = 1
        self.n_baryon = 8
        self.C.phi_meson = 0
        self.C.hyper_phi_sigma = 1
        self.filenames['mass_crust'] = 'mass_hyper.dat'
        self.filenames['mass_nocrust'] = None
        self.filenames['eos'] = 'hyper.dat'
        self.filenames.update(etap_f='etap_f.dat')
        self.filenames.update(etap_n='etap_n.dat')
        self.filenames['meff'] = 'meff_hyper.dat'
        self.filenames['vs'] = 'vs_hyper.dat'

    def dumpEos(self):
        E, P, n = self.EPN()
        tab = [self.nrange/self.n0, E, P]
        conc = self.concentrations()
        tab += [conc[:,i] for i in range(self.n_baryon)]
        tab += [self.rho[:,0]]
        table = tabulate(arr(tab).transpose(), ['n/n0', 'E_NS[m_\pi^4]',
                                                'P_NS[m_\pi^4]'] +
        self.part_names[:self.n_baryon] + ['f'], tablefmt='plain'
        )
        with open(join(self.foldername, self.filenames['eos']), 'w') as f:
            f.write(table)

        return tab

    def dumpEtap(self):
        frange = np.linspace(0, 1., 100)
        tab = arr([frange, map(lambda z: self.C.eta_p(z), frange)]).transpose()
        with open(join(self.foldername, self.filenames['etap_f']), 'w') as f:
            f.write(tabulate(tab, ['f', 'eta_p(f)'], tablefmt='plain'))

    def dumpChi(self):
        pass


class HyperonPhi(Hyperon):
    def __init__(self, C):
        super(HyperonPhi, self).__init__(C)
        self.C.Hyper = 1
        self.C.phi_meson = 1
        self.C.phi_kind = 1
        self.C.hyper_sigma_kind = 1
        self.filenames['mass_crust'] = 'mass_hyper_phi.dat'
        self.filenames['mass_nocrust'] = None
        self.filenames['eos'] = 'hyper_phi.dat'
        self.filenames['etap_f'] = 'etap_phi_f.dat'
        self.filenames.update(etap_n='etap_phi_n.dat')
        self.filenames['meff'] = 'meff_hyper_phi.dat'
        self.filenames['vs'] = 'vs_hyper_phi.dat'






class HyperonPhiSigma(HyperonPhi):
    def __init__(self, C):
        super(HyperonPhiSigma, self).__init__(C)
        self.C.Hyper = 1
        self.C.phi_meson = 1
        self.C.sigma_kind = 1
        self.C.hyper_sigma_kind = 0
        self.filenames['mass_crust'] = 'mass_hyper_phi_sigma.dat'
        self.filenames['mass_nocrust'] = None
        self.filenames['eos'] = 'hyper_phi_sigma.dat'
        self.filenames['etap_f'] = 'etap_phi_sigma_f.dat'
        self.filenames.update(etap_n='etap_phi_sigma_n.dat')
        self.filenames['meff'] = 'meff_hyper_phi_sigma.dat'
        self.filenames['vs'] = 'vs_hyper_phi_sigma.dat'

class Model(Wrapper):
    def __init__(self, C, K0=None, f0=None, J0=None, suffix=None):
        if any([K0, f0, J0]):
            print map(eval, [name for name in ['K0', 'f0', 'J0']])
            names = ['K0', 'f0', 'J0']
            actual = [name for name in names if eval(name) is not None]
            print actual
            args = dict(zip(actual, map(eval, actual)))
            print args
            wr = Wrapper_old(C())
            wr.solve(**args)
            params = ['Cs', 'Co', 'Cr', 'b', 'c', 'f0']
            pvals = [eval('wr.C.' + item) for item in params]
        super(Model, self).__init__(C)
        self.nucl = Nucleon(C)
        self.sym = Sym(C)
        self.neutr = Neutr(C)
        self.hyper = Hyperon(C)
        self.hyper_phi = HyperonPhi(C)
        self.hyper_phi_sigma = HyperonPhiSigma(C)
        self.children = [self.nucl, self.sym, self.neutr, self.hyper,
                         self.hyper_phi, self.hyper_phi_sigma]

        if any([K0, f0, J0]):
            for subwr in [self] + self.children:
                [six.exec_(s) for s in ['subwr.C.' + item + ' = pvals[' + str(i) + ']'
                           for i, item in enumerate(params)]]

        if suffix is not None:
            self.foldername = os.path.join(BASEFOLDER,
                                           self.Ctype.__name__.strip('_') +
                                           suffix)
            if not os.path.exists(self.foldername):
                os.makedirs(self.foldername)
            for child in self.children:
                child.foldername = self.foldername
        self.filenames['eos'] = 'eos.dat'
        self.filenames.update(eta='eta_F.dat')
        self.filenames.update(props='props.dat')
        self.filenames.update(intro='intro.dat')
        self.mpi3tofmm3 = (self.m_pi/197.33)**3







    def dumpEos(self, nmax=None, npoints=None, write=True):
        if nmax is None:
            nmax = self.nmax
        if npoints is None:
            npoints = self.npoints

        nrange = self.nrange

        esym, fsym = self.sym.Ebind(nrange, ret_f=True)
        psym = self.sym.P(nrange)

        en, fn = self.neutr.Ebind(nrange, ret_f=True)
        pn = self.neutr.P(nrange)

        E, P, n = self.nucl.EPN()
        conc = self.nucl.concentrations()
        n_e, n_mu, x_du = self.nucl.lepton_concentrations(ret_du=1)
        table = arr([nrange/self.n0, esym, psym, fsym,
                     E, P/self.mpi4_2_mevfm3, conc[:, 1], self.nucl.rho[:, 0],
                     en, pn, fn, self.nucl.mu_e*self.m_pi, x_du]).transpose()

        tab = tabulate(table, ['n/n_0', 'e_{symm} [MeV]',
                               'P_{symm} [MeV / fm^3] ',
                               'f(n){symm}',
                               'E_NS [m_\pi^4]',
                               'P_NS, [m_\pi^4]',
                               'n_p/n', 'f(n){NS}',
                               'E_N [MeV]', 'P_N [MeV/fm^3]',
                               'f(n){Pure NM}',
                               'mu_e [MeV]',
                               'x_DU'],
                       tablefmt='plain')

        # print tab
        if write:
            with open(join(self.foldername, self.filenames['eos']), 'w') as f:
                f.write(tab)

        f = open(join(self.foldername, 'for_grigorian.dat'), 'w')
        _n = self.nucl.nrange * self.mpi3tofmm3
        mu_n = map(lambda z: self.m_pi*eos.mu(z, 1, self.C),
                   self.nucl.rho)
        meff = map(lambda z: self.C.phi_n(0, z), self.nucl.rho[:, 0])
        table = arr([mu_n, P, E*self.mpi4_2_mevfm3,
                     _n,
                     _n * conc[:, 0],
                     _n * conc[:, 1],
                     meff,
                     meff,
                     n_e,
                     n_mu,
                     np.zeros(_n.shape),
                     np.zeros(_n.shape),
                     np.zeros(_n.shape),
                     self.nucl.mu_e * self.m_pi,
                     conc[:, 1]
                    ]).transpose()

        tab = tabulate(table, ['mu_n [MeV]', 'P [MeV/fm^3]', 'E[MeV/fm^3]',
                               'n_b [1/fm^3]',
                               'n_n [1/fm^3]', 'n_p [1/fm^3]', 'm_n^*', 'm_p^*',
                               'n_e [1/fm^3]', 'n_mu [1/fm^3]', 'n_u [1/fm^3]',
                               'n_d [1/fm^3]', 'n_s [1/fm^3]',
                               'mu_e [MeV]', 'Y_p'], floatfmt='.6f')
        f.write(tab)
        f.close()
        return table

    def J(self, nrange):
        return self.neutr.Ebind(nrange) - self.sym.Ebind(nrange)

    def Jtilde(self, nrange):
        return self.sym.Jtilde(nrange=nrange)

    def dumpAll(self, hyper=True):
        self.dumpEos()
        self.nucl.dumpMeff()
        self.nucl.dumpMassesCrust()
        self.nucl.dumpEos()
        for s in [self.nucl, self.sym, self.neutr]:
            s.dumpVs()
            s.dumpParts()

        mods = []
        if hyper:
            mods = [self.hyper, self.hyper_phi, self.hyper_phi_sigma]
            # mods = [self.hyper_phi_sigma]
        for s in mods:
            s.dumpEos()
            s.dumpMassesCrust()
            s.dumpMeff()
            s.dumpEtap()

    def dumpInspect(self):
        pass

    def dumpScalings(self):
        tabF = []
        C = self.C
        for f in np.linspace(0., 1., 100., endpoint=False):
            tabF.append([f, C.eta_s(f) + 2*C.Cs*C.U(f)/(C.M[0]**4 * f**2),
                          C.eta_o(f), C.eta_r(f),
                          C.eta_p(f), C.phi_n(0, f), C.U(f), C.eta_s(f)])
        tableF = tabulate(tabF, headers=['f',
                                 'eta_sigma',
                                 'eta_omega',
                                 'eta_rho',
                                 'eta_phi',
                                 'Phi_N',
                                 'U',
                                 'eta_sigma_bare'],tablefmt='plain', floatfmt='.10f')
        with open(join(self.foldername, self.filenames['eta']), 'w') as f:
            f.write(tableF)

    def dumpProps(self):
        tabParams= [['Cs', self.C.Cs],
         ['Co', self.C.Co],
         ['Cr', self.C.Cr],
         ['b', self.C.b],
         ['c', self.C.c],
         ['Csp', self.C.Csp]
         ]

        tabPart = [['m_N', self.m_pi*self.C.M[0]],
                   ['m_Lambda', self.m_pi*self.C.M[2]],
                   ['m_Sigma', self.m_pi*self.C.M[4]],
                   ['m_Xi', self.m_pi*self.C.M[6]]
                   ]

        tabRes = [['E_bind [MeV]',eos.EBind(np.array([self.C.f0, self.n0/2,
                                                      self.n0/2]), self.C)],
                  ['K [MeV]', self.K()],
                  ["K' [MeV]", self.Kprime()],
                  ['J [MeV]', self.J()],
                  ['L [MeV]', self.L()],
                  ['Ksym [MeV]', self.Ksymm()],
        ]

        Params = tabulate(tabParams,floatfmt='.10f',tablefmt='plain')
        Part = tabulate(tabPart,tablefmt='plain',floatfmt='.10f')
        Res = tabulate(tabRes,tablefmt='plain',floatfmt='.10f')
        Xs = tabulate([[self.part_names[i], self.C.X_s[i]] for i in range(8)], tablefmt='plain')

#         cPar = ''
#         for member in inspect.getmembers(self.C, inspect.isdatadescriptor):
#             print member
#
#         exit()
        f = open(os.path.join(self.foldername, self.filenames['props']), 'w')
        f.write(Params + '\n' + Part + '\n' + Res + '\n' + 'X_\sigma\n' + Xs + '\n')
        f.close()

    def K(self):
        return eos.K(self.n0, self.C)

    def J(self):
        return eos.J(self.n0, self.C)

    def L(self):
        return 3*self.n0*derivative(lambda z: eos.J(z, self.C),
                                    self.n0, dx=1e-3)

    def Kprime(self):
        return -3*self.n0*(derivative(lambda z: eos.K(z, self.C),
                                      self.n0, dx=1e-3, order=3) -
                2*eos.K(self.n0,self.C)/self.n0)

    def Ksymm(self):
        return 9*self.n0**2 * derivative(lambda z: eos.J(z, self.C),
                                         self.n0, dx=1e-3, n=2)

    def UofE(self, i, n):
        f = eos.f_eq(n, np.array([self.C.f0]), 1, self.C)
        pots = eos.potentials(np.insert(n, 0, f), 5, self.C)
        print pots
        V = self.C.X_o[i] * pots[2]
        S =  (-self.C.M[i]) * (1 - self.C.phi_n(i, self.C.X_s[i]*self.C.M[0]/self.C.M[i]*f[0]))
        Uopt = lambda z: z + self.C.M[0] - sqrt((z + self.C.M[0] - V)**2 - S*(2 * self.C.M[i] + S))
        erange = np.linspace(-50./self.m_pi, 10., 100)
        return self.m_pi*erange, self.m_pi*Uopt(erange)

    def dumpUn(self, nmin=0., nmax=5., npoints=100):
        nrange = np.linspace(nmin, nmax, npoints)
        res = []
        i = 0
        for n in nrange:
            n_in = np.array([n/2, n/2])
            f = eos.f_eq(n_in, np.array([self.C.f0]), 1, self.C)
            pots = eos.potentials(np.insert(n_in, 0, f), 5, self.C)
            print pots
            V = self.C.X_o[i] * pots[2]
            S =  (-self.C.M[i]) * (1 - self.C.phi_n(i, self.C.X_s[i]*self.C.M[0]/self.C.M[i]*f[0]))
            Uopt = self.m_pi*(V + S)
            res.append(Uopt)
#         plt.plot(nrange/self.n0, res)
#         plt.show()
        tab = np.array([nrange/self.n0, np.array(res)]).transpose()
        table = tabulate(tab, ['n/n_0', 'U_N [MeV]'], tablefmt='plain')
        with open(join(self.foldername, 'U_N.dat'), 'w') as f:
            f.write(table)
        return nrange, np.array(res)

    def dumpUofE(self, show=False):
        nU = np.array([self.n0/2, self.n0/2])
        if show:
            lines = []
        ulist = []
        for i in [0, 2, 3, 6]:
            e, u = self.UofE(i, nU)
            if i == 0:
                ulist.append(e)
            ulist.append(u)
            if show:
                line, = plt.plot(e, u)
                lines.append(line)
        if show:
            plt.legend(lines, ['N', 'L', 'S', 'X'])
            plt.show()

        ulist = np.array(ulist).transpose()
        tab = tabulate(ulist, ['E [MeV]', 'U_N [MeV]', 'U_L [MeV]',
                               'U_S, [Mev]', 'U_X [MeV]'], tablefmt='plain')

        with open(join(self.foldername, 'UofE.dat'), 'w') as f:
            f.write(tab)

    def dumpUofE_anti(self, show=False):
        nU = np.array([self.n0/2, self.n0/2])
        if show:
            lines = []
        ulist = []
        for i in [0, 2, 3, 6]:
            e, u = self.UofE_anti(i, nU)
            if i == 0:
                ulist.append(e)
            ulist.append(u)
            if show:
                line, = plt.plot(e, u)
                lines.append(line)
        if show:
            plt.legend(lines, ['N', 'L', 'S', 'X'])
            plt.show()

        ulist = np.array(ulist).transpose()
        tab = tabulate(ulist, ['E [MeV]', 'U_N [MeV]', 'U_L [MeV]',
                               'U_S, [Mev]', 'U_X [MeV]'], tablefmt='plain')

        with open(join(self.foldername, 'UofE_anti.dat'), 'w') as f:
            f.write(tab)

    def UofE_anti(self, i, n):
        f = eos.f_eq(n, np.array([self.C.f0]), 1, self.C)
        pots = eos.potentials(np.insert(n, 0, f), 5, self.C)
        print pots
        V = -self.C.X_o[i] * pots[2]
        S =  (-self.C.M[i]) * (1 - self.C.phi_n(i, self.C.X_s[i]*self.C.M[0]/self.C.M[i]*f[0]))
        Uopt = lambda z: z + self.C.M[0] - sqrt((z + self.C.M[0] - V)**2 - S*(2 * self.C.M[i] + S))
        erange = np.linspace(-50./self.m_pi, 10., 100)
#         plt.plot(erange, self.m_pi*Uopt(erange))
#         plt.show()
        return self.m_pi*erange, self.m_pi*Uopt(erange)

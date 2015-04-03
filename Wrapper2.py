import eosWrap as eos
import numpy as np
from numpy import array as arr, pi, sqrt
import matplotlib.pyplot as plt
import os
from tabulate import tabulate
from scipy import interpolate

BASEFOLDER = '/home/const/wr2_test'


class Output(object):
    def __init__(self, foldername):
        self.foldername = foldername

    def plot(self, f):
        pass

    def write(self, f, filename):
        pass


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
        self.nmin = 0.
        self.nmax = 4.
        self.n_baryon = 2
        self.verbose = False
        self.mpi4_2_mevfm3 = self.m_pi * (self.m_pi / 197.33) ** 3
        self.nmin = 0.
        self.nmax = 10*self.n0
        self.npoints = 1000
        self.nrange = np.linspace(self.nmin, self.nmax, self.npoints,
                                  endpoint=False)



    def reset(self, iterations=30, timeout=None):
        """Calculates the EoS stored in the wrapper. Has to be overrided in
        Symm, Neutron. 
        
        Calculates the particle composition, scalar field, energy density, 
        pressure, ... for the given parameter set self.C. Wrapper.set is 
        set to True after s succesful calculation. 
        
        Note: \sigma^* is not supported in the current version.

        """
        init = arr([0. for i in range(self.n_baryon - 1)])

        rho = []
        mu_e = []
        f = arr([0.])  # sp omitted
        for i, _n in enumerate(self.nrange):
            if timeout is None:
                init = eos.stepE(_n, init, f, len(init), iterations, self.C)
            else:
                pass

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
        self._E = np.array(map(lambda z: eos.E(z, self.C), self.rho))
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
                if np.allclose(nrange, self.nrange):
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


    def dumpMassesCrust(self, nmin=0.6, nmax=4., npoints=100, write=False):
        pass

    def stars_crust(self, ncut_crust=0.45, ncut_eos=0.6, inter='cubic',
                    n_stars=None, nmin=.6, nmax=5.0, npoints=50,
                    crust="crust.dat", show=False, crustName=None,
                    ret_str=False, fasthyp=False, neutron=0):
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

        #interpolating hyperon fractions
        if self.C.Hyper:
            hyp_mult = np.array([0, 0, 1, 1, 1, 1, 2, 2])
            inter_hyp = [interpolate.interp1d(self.n, hyp_mult[i] *
                                              self.concentrations()[:, i])
                         for i in range(2, 8)]

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
            if self.C.Hyper and ret_str:
                hyper_N = []
                hyper_NR = []
                for f in inter_hyp:
                    hyper_N.append(lastN * map(f, lastN))
                    hyper_NR.append(lastN * map(f, lastN) * grav_mult)
                #                 if lastM[-1] > 1.6:
                #                     plt.plot(lastN/self.n0, np.array(hyper_N).transpose())
                #                     plt.show()
                res_str = 0
                for _N in hyper_N:
                    res_str += np.multiply(_N, grav_mult)
                #                 if lastM[-1] > 1.6:
                #                     plt.plot(lastR, grav_mult)
                #                     for f in inter_hyp:
                #                         plt.plot(lastR, grav_mult * np.array(map(f, lastN)))
                #                     plt.show()

                str_frac = np.trapz(res_str, dx=dx) / (3 * integral)

                str_fracs.append(str_frac)

            Mgrav.append((0.0004898007281478712) * integral)

        MR = np.array(MR)
        Mgrav = np.array(Mgrav)
        str_fracs = np.array(str_fracs)

        if ret_str:
            return (n_stars, MR[:, 0], MR[:, 1], 931.5 / self.m_pi * MR[:, 2],
                    931.5 / self.m_pi * Mgrav, str_fracs)
        else:
            return (n_stars, MR[:, 0], MR[:, 1], 931.5 / self.m_pi * MR[:, 2],
                    931.5 / self.m_pi * Mgrav)


    def dumpScalings(self):
        pass

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
        return np.array(res)


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
        for n in nlist:
            f = eos.f_eq(arr(n), arr([f]), 1, self.C)[0]
            flist.append(f)
            elist.append(eos._E(np.insert(n, 0, f), self.C))
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


class Nucleon(Wrapper):
    def __init__(self, C):
        super(Nucleon, self).__init__(C)
        self.C.Hyper = 0

    def E(self, nrange, ret_f=False, f=0.):
        E, P, n = self.EPN(nrange=nrange)
        if ret_f:
            flist = self.rho[:, 0]
            return E, flist
        else:
            return E


    def lepton_concentrations(self, ret_du=False):
        """Returns lepton concentrations for self.nrange. """
        self.check()
        ne_list = []
        nmu_list = []
        du_list=[]
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





class Sym(Nucleon):
    def __init__(self, C):
        super(Sym, self).__init__(C)
        self.C.Hyper = 0

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


class Neutr(Nucleon):
    def __init__(self, C):
        super(Neutr, self).__init__(C)
        self.C.Hyper = 0

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


class HyperonPhi(Hyperon):
    def __init__(self, C):
        super(HyperonPhi, self).__init__(C)
        self.C.Hyper = 1
        self.C.phi_meson = 1


class HyperonPhiSigma(HyperonPhi):
    def __init__(self, C):
        super(HyperonPhi, self).__init__(C)
        self.C.Hyper = 1
        self.C.phi_meson = 1


class Model(Wrapper):
    def __init__(self, C):
        super(Model, self).__init__(C)
        self.nucl = Nucleon(C)
        self.sym = Sym(C)
        self.neutr = Neutr(C)
        self.hyper = Hyperon(C)
        self.hyper_phi = HyperonPhi(C)



    def dumpEos(self, nmax=None, npoints=None):
        if nmax is None:
            nmax = self.nmax
        if npoints is None:
            npoints = self.npoints

        nrange = self.nrange

        Esym, fsym = self.sym.Ebind(nrange, ret_f=True)
        Psym = self.sym.P(nrange)

        En, fn = self.neutr.Ebind(nrange, ret_f=True)
        Pn = self.neutr.P(nrange)

        E, P, n = self.nucl.EPN()
        conc = self.nucl.concentrations()
        n_e, n_l, x_du = self.nucl.lepton_concentrations(ret_du=1)
        table = arr([nrange/self.n0, Esym, Psym, fsym,
                     E, P, conc[:,1], self.nucl.rho[:,0],
                     En, Pn, fn, self.nucl.mu_e*self.m_pi, x_du]).transpose()

        tab = tabulate(table, ['n/n_0', 'e_{symm} [MeV]',
                               'P_{symm} [MeV / fm^3] ',
                               'f(n){symm}' ,
                               'E_NS [m_\pi^4]',
                               'P_NS, [m_\pi^4]',
                               'n_p/n', 'f(n){NS}',
                               'E_N [MeV]','P_N [MeV/fm^3]',
                               'f(n){Pure NM}',
                               'mu_e [MeV]',
                               'x_DU'],
                       tablefmt='plain')

        print tab

















        

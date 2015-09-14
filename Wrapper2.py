import eosWrap as eos
# from lmfit.minimizer import minimize
# from lmfit.parameter import Parameters
from matplotlib.widgets import Slider
import numpy as np
from numpy import array as arr, pi, sqrt
import matplotlib.pyplot as plt
import os
from os.path import join
from scipy.interpolate.interpolate import interp1d
from scipy.optimize.minpack import leastsq
from tabulate import tabulate
from scipy import interpolate, optimize
from scipy.misc import derivative
from multiprocessing import Queue, Process
from Wrapper import Wrapper as Wrapper_old
import six


BASEFOLDER = '/home/const/MKVdelta_local/data_new/'

class Wrapper(object):
    def __init__(self, Ctype, basefolder_suffix=''):
        self.Ctype = Ctype
        self.C = Ctype()
        self.foldername = os.path.join(BASEFOLDER + basefolder_suffix, Ctype.__name__.strip('_'))
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
        self.mpi3tofmm3 = (self.m_pi/197.33)**3
        self.nmin = 0.
        self.nmax = 10*self.n0
        self.npoints = 400
        self.nrange = np.linspace(self.nmin, self.nmax, self.npoints,
                                  endpoint=False)
        #example of setting filenames for a descendant of Wrapper
        self.filenames = {'mass_nocrust': None, 'eos': None,
                          'mass_crust' : None, 'meff':None,
                          'eta' : None, 'vs' : None,
                          'Eparts':None,
                          'Pparts':None,
                          'uniparts': None,
                          'profiles' : None,
                          'grig' : None,
                          'mu' : None}

        self.part_names = ['n', 'p', 'Lambda', 'Sigma-',
                           'Sigma0', 'Sigma+', 'Xi-', 'Xi0', 'e', 'mu']

        self.md = MassDriver()


    def stepE(self, n, last, f, lenLast, iter, C, que):
        rho = eos.stepE(n, last, f, len(last), iter, C)
        que.put(rho)

    def reset(self, iterations=100, timeout=20):
        """Calculates the EoS stored in the wrapper. Has to be overriden in
        Symm, Neutron. 
        
        Calculates the particle composition, scalar field, energy density, 
        pressure, ... for the given parameter set self.C. Wrapper.set is 
        set to True after s succesful calculation. 
        
        Note: \sigma^* is not supported in the current version.

        """

        print("Resetting " + self.__repr__() + " for model " + self.Ctype.__name__.strip('_'))
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
                    print("timeout reached")
                    self.rho = np.ascontiguousarray(rho[:])
                    self.mu_e = np.ascontiguousarray(arr(mu_e))
                    self.nrange = self.nrange[: self.rho.shape[0]]
                    _E = [eos.E(z, self.C) for z in self.rho]
                    self._E = np.ascontiguousarray(np.array(_E[:]))
                    self._P = np.ascontiguousarray(self.P_chem(self.rho))
                    self.set=1
                    return
                init = queue.get(timeout=None)
                # print init

            if i % (len(self.nrange) / 20) == 0:
                print('.', end=' ')

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
        dEparts = []
        for part in self.Eparts.transpose():
            dEparts.append(np.gradient(part, self.nrange[1]-self.nrange[0]))
        dEparts = arr(dEparts)
        p1 = self.nrange * dEparts
        self.Pparts = p1.transpose() - self.Eparts
        # self._E = np.array(map(lambda z: eos.E(z, self.C), self.rho))
        self._E = np.ascontiguousarray(self._E)
        self._P = np.ascontiguousarray(self.P_chem(self.rho))
        self.set = True

    def EPN(self, nrange=None):
        if self.check(nrange=nrange):
            return np.copy(self._E), np.copy(self._P), np.copy(self.nrange)
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

    # def dumpMasses(self, nmin=0.4, nmax=4., npoints=100, write=False):
    #     if write:
    #         if not os.path.exists(self.foldername):
    #             os.makedirs(self.foldername)
    #     NPOINTS_EOS = 1000
    #     if self.check(nmin=0., nmax=nmax, npoints=NPOINTS_EOS):
    #         if not self.driverSet:
    #             self.setDriver()
    #
    #     n, m, r, mb = self.stars(npoints=100, nmax=10. * self.n0)
    #     if write:
    #         table = np.array([n / self.n0, m, r, mb]).transpose()
    #         f = open(os.path.join(self.foldername, 'masses.dat'), 'w')
    #         tab = tabulate(table, ['n/n_0',
    #                                'M [M_{sun}]',
    #                                'R [km]',
    #                                'M_B [M_sun]'],
    #                        tablefmt='plain')
    #         f.write(tab)
    #         f.close()
    #     return n, m, r


    def dumpMassesCrust(self, nmin=0.4, nmax=4.8, npoints=100, write=True, fname=None, ret_frac=True):
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

    def dumpMasses(self, nmin=.4, nmax=5., npoints=100, write=True, fname=None, ret_frac=True):
        self.check()
        out = self.stars_nocrust(nmin=nmin, nmax=nmax, npoints=npoints, ret_frac=ret_frac)
        if write:
            tab = arr(out).transpose()
            names = ['n/n0', 'M[M_sun]', 'R[km]']
            if ret_frac:
                names += self.part_names[1:self.n_baryon]
            table = tabulate(tab, names, tablefmt='plain')
            if fname is None:
                fname = self.filenames['mass_nocrust']
            with open(join(self.foldername, fname), 'w') as f:
                f.write(table)
        return out

    def stars_nocrust(self, nmin=.6, nmax=5., npoints=50, show=False, ret_frac=False):
        E, P, N = self.EPN(self.nrange)
        P = P/self.mpi4_2_mevfm3
        # E *= self.m_pi ** 4
        # P *= self.m_pi ** 4
        self.dr = eos.KVDriver()
        self.dr.set(E, P, N)
        MR = []
        nstars = np.linspace(nmin, nmax, npoints)
        for _n in nstars:
            MR.append(eos.star_crust2(_n, 3, self.dr, 1e-11))
        MR = np.array(MR)
        return [nstars, MR[:, 0], MR[:, 1]]

    def stars_crust(self, ncut_crust=0.45, ncut_eos=0.7, inter='cubic',
                    n_stars=None, nmin=.6, nmax=5.0, npoints=50,
                    crust="crust.dat", show=False, crustName=None,
                    ret_frac=False, fasthyp=False, neutron=0, ret_i=0, force_reset=0):
        neos = self.npoints
        if n_stars is not None:
            nmax = n_stars[-1]  # TODO too dirty workaround
        if not np.all([self.md.E, self.md.N, self.md.P]) or force_reset:
            self.check(nrange=self.nrange)
            E, P, N = self.EPN(self.nrange)
            P = P / self.mpi4_2_mevfm3
            N = N / self.n0
            print(self.Ctype.__name__)
            np.savetxt('eos'+self.Ctype.__name__+'.dat', arr([E, P, N]).transpose())
            E *= self.m_pi ** 4
            P *= self.m_pi ** 4

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

            crust = np.loadtxt("/home/const/workspace/swigEosWrapper/crust.dat")
            crust[:, 0] /= self.m_pi**4
            crust[:, 1] /= self.m_pi**4
            np.savetxt( 'crust_export.dat', crust)
            n_eos = 5
            i_n_eos = np.argmin(abs(N - [ncut_eos for i in N]))



            plist = np.append(p[:], P[i_n_eos:(i_n_eos + n_eos)])
            elist = np.append(e[:], E[i_n_eos:(i_n_eos + n_eos)])
            nlist = np.append(n[:], N[i_n_eos:(i_n_eos + n_eos)])

            print(i_n_eos)
            print(P[i_n_eos:(i_n_eos + n_eos)])
            print(N[i_n_eos:(i_n_eos + n_eos)])
            # exit()
            iP = interpolate.interp1d(nlist, plist, kind=inter)
            iE = interpolate.interp1d(nlist, elist, kind=inter)

            gamma = 1. / 4.
            iN = np.linspace(0, ncut_eos ** gamma, 10000)
            iN = iN ** (1. / gamma)

            crust_p = np.array(list(map(iP, iN)))
            crust_e = np.array(list(map(iE, iN)))

            #         finalE = np.append(crust_e, E[i_n_eos+n_eos:])
            #         finalP = np.append(crust_p, P[i_n_eos + n_eos:])
            #         finalN = np.append(iN, N[i_n_eos+n_eos:])

            finalE = np.append(crust_e, E[i_n_eos + n_eos:]) / self.m_pi ** 4
            finalP = np.append(crust_p, P[i_n_eos + n_eos:]) / self.m_pi ** 4
            finalN = np.append(iN, N[i_n_eos + n_eos:])
            finalE[0] = 0
            finalP[0] = 0
            print(finalN)
            self.md.setEos(finalN, finalE, finalP)
        else:
            finalE = self.md.E
            finalN = self.md.N
            finalP = self.md.P

        # plt.plot(n, arr(p)/self.m_pi**4, label='crust')
        # plt.plot(crust[:,2], arr(crust[:, 1]), label='crust_whole')
        # plt.plot(N, arr(P)/self.m_pi**4, label='eos')
        # plt.plot(finalN, finalP, label='interp')
        # plt.legend()
        # plt.xlabel('n/n_0')
        # plt.ylabel('P [m_\pi^4]')
        # plt.xlim([0, 1])
        # plt.ylim([0., 0.05])
        # plt.show()

        mpi2km = 5.7202e-5
        np.savetxt('EPkm2_KVOR06.dat',
                   np.array([finalN, mpi2km*finalE, mpi2km*finalP]).transpose(), fmt='%.6e', header='n/n0  E  P')

        np.savetxt('EPmpi4_KVOR06.dat',
                   np.array([finalN, finalE, finalP]).transpose(), fmt='%.6e', header='n/n0  E  P')

        # plt.plot(finalE, finalE)
        # plt.show()
        #
        # plt.plot(np.diff(finalE))
        # plt.show()

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

        self.dr = eos.KVDriver(finalE, finalP, finalN*self.n0)
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
            print(_n)
            MR.append(eos.star_crust2(_n, 3, self.dr, 1e-10))
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
                    hyper_N.append(lastN * list(map(f, lastN)))
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
            return [n_stars, MR[:, 0], MR[:, 1], 931.5 / self.m_pi * MR[:, 2],
                    931.5 / self.m_pi * Mgrav]


    def dumpScalings(self):
        E, f = self.E(self.nrange, ret_f=1)
        C = self.C
        tab = arr([[z, C.eta_s(z), C.eta_o(z), C.eta_r(z), C.phi_n(0, z), C.U(z)] for z in f])
        table = tabulate(np.insert(tab, 0, self.nrange/self.n0, axis=1), ['n/n0', 'f', 'eta_s', 'eta_o',
                                                                          'eta_r', 'phi_n', 'U'],
                         tablefmt='plain')
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


    def E_gen(self, nlist, ret_f=False, f=0., dn=None):
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
        dEparts = []
        for n in nlist:
            epart = np.zeros(9, dtype='float')
            f = eos.f_eq(arr(n), arr([f]), 1, self.C)[0]
            flist.append(f)
            elist.append(eos._E(np.insert(n, 0, f), self.C, epart))
            eparts.append(epart)
            self.Eparts = arr(eparts)
            if dn is not None:

                # dEparts.append(np.gradient(part, sum(nlist[1]) - sum(nlist[0])))
                de1 = np.zeros(9, dtype='float')
                de2 = np.zeros(9, dtype='float')
                eos._E(np.insert(n+dn/2, 0, f), self.C, de1)
                eos._E(np.insert(n-dn/2, 0, f), self.C, de2)
                dEparts.append((arr(de1) - arr(de2))/np.sum(dn, axis=0))
        if dn is not None:
            dEparts = arr(dEparts)

            p1 = np.sum(nlist, axis=1) * dEparts.transpose()
            self.Pparts = p1.transpose() - self.Eparts
        if ret_f:
            return [arr(elist), arr(flist)]
        else:
            return arr(elist)


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
            print('Concentrations: wrapper is not set! Pass.')
            self.reset()
            return

        rho = []
        sp = 1 + self.C.sprime
        for _r in self.rho:
            if np.sum(_r[sp:]) > 0:
                rho.append(_r[sp:] / np.sum(_r[sp:]))
            else:
                rho.append([1.] + [0. for i in _r[sp + 1:]])

        return np.array(rho)

    def dumpVs(self, nrange=None, write=True):
        E, P, n = self.EPN(nrange=nrange)
        vs = np.gradient(P) / np.gradient(E) / self.mpi4_2_mevfm3
        if write:
            with open(join(self.foldername, self.filenames['vs']), 'w') as f:
                f.write(tabulate(arr([n/self.n0, vs]).transpose(), ['n/n0', 'v_s^2'], tablefmt='plain'))
        return vs

    def dumpParts(self):
        self.check()
        if self.filenames['Eparts']:
            table = tabulate(np.insert(self.Eparts*self.mpi4_2_mevfm3, 0, self.nrange/self.n0, axis=1),
                             ['n/n0', 'f', 'U(f)', 'K_n', 'K_p', 'om', 'rho',
                              'phi', 'e', 'mu'], tablefmt='plain')
            with open(join(self.foldername, self.filenames['Eparts']), 'w') as f:
                f.write(table)

        if self.filenames['Pparts']:
            table = tabulate(np.insert(self.Pparts*self.mpi4_2_mevfm3, 0, self.nrange/self.n0, axis=1),
                             ['n/n0', 'f', 'U(f)', 'K_n', 'K_p', 'om', 'rho',
                              'phi', 'e', 'mu'], tablefmt='plain')
            with open(join(self.foldername, self.filenames['Pparts']), 'w') as f:
                f.write(table)

        if self.filenames['uniparts']:
            E, P, n = self.EPN()
            Ep = self.Eparts * self.mpi4_2_mevfm3
            Pp = self.Pparts * self.mpi4_2_mevfm3
            eparts = [n/self.n0,
                          E * self.mpi4_2_mevfm3,
                          Ep[:, 0] + Ep[:, 1],
                          Ep[:, 2], Ep[:, 3],
                          Ep[:, 4] + Ep[:, 5]]

            pparts =     [P,
                          Pp[:, 0] + Pp[:, 1],
                          Pp[:, 2], Pp[:, 3],
                          Pp[:, 4] + Pp[:, 5]]

            table = tabulate(arr(eparts + pparts).transpose(),
                             ['n/n0', 'E', 'Esigma', 'Ekn', 'Ekp', 'Evector',
                              'P', 'Psigma', 'Pkn', 'Pkp', 'Pvector'], tablefmt='plain')
            with open(join(self.foldername, self.filenames['uniparts']), 'w') as f:
                f.write(table)

    def getProfile(self, n_c, ret_frac=False):
        self.dumpMassesCrust(nmin=n_c, nmax=n_c, npoints=1, write=False)
        lastN = self.dr.getLastN(self.dr.nSize)[:-1]
        lastR = self.dr.getLastR(self.dr.nSize)[:-1]
        lastM = self.dr.getLastM(self.dr.nSize)[:-1]
        grav_mult = []
        for i, r in enumerate(lastR):
            grav_mult.append(r**2. / np.sqrt(1 - 2 * 1.4677 * lastM[i] / r))
        ###This grav_mult differs by a lack of r^2 from that in dumpMassesCrust
        ###Already not
        grav_mult = np.array(grav_mult)
        dx = lastR[1] - lastR[0]
        mb = []

        for i, _r in enumerate(lastR):
            mb.append(931.5 / self.m_pi * 0.0004898007281478712*
                      np.trapz(np.multiply(lastN[:i], np.array(grav_mult[:i])), dx=dx))

        conc = self.concentrations()

        if not ret_frac:
            return lastN, lastR, lastM, grav_mult #, arr([[f(z) for f in inter_hyp] for z in lastN]).transpose().tolist()
        else:
            inter_hyp = [interpolate.interp1d(self.nrange, conc[:, i])
            for i in range(1, self.n_baryon)]
            return lastN, lastR, lastM, np.array(mb), \
                   arr([[f(z) for f in inter_hyp] for z in lastN]).transpose().tolist()

    def dumpProfiles(self, n_c):
        n, r, mg, gm, part = self.getProfile(n_c, ret_frac=1)
        # n, r, mg, gm= self.getProfile(n_c)
        E, P, n_eos = self.EPN()
        P /= self.mpi4_2_mevfm3
        vs = self.dumpVs(write=0)
        gamma = np.nan_to_num((1 + E / P) * vs)
        iGamma = interp1d(n_eos, gamma)
        with open(join(self.foldername,    '%.3f'%(n_c/self.n0)+self.filenames['profiles']), 'w') as f:
            f.write(tabulate(arr([n/self.n0, r] + part + [mg] + [gm] + [iGamma(n)]).transpose(), ['n/n0', 'r[km]']
                           + self.part_names[1:] + ['M', 'Mb', 'Gamma'], tablefmt='plain'))
            # f.write(tabulate(arr([n/self.n0, r] + [mg] + [gm] + [iGamma(n)]).transpose(), ['n/n0', 'r[km]', 'M[M_sun]']
            #                + ['gamma', 'Gamma'], tablefmt='plain', floatfmt='.3f'))


    def dumpGrig(self, E, P, conc, n_e, n_mu):
        f = open(join(self.foldername, self.filenames['grig']), 'w')
        _n = self.nrange * self.mpi3tofmm3
        mu_n = [self.m_pi * eos.mu(z, 1, self.C) for z in self.rho]
        meff = [self.C.phi_n(0, z) for z in self.rho[:, 0]]
        table = arr([mu_n, P, E * self.mpi4_2_mevfm3,
                     _n,
                     _n * conc[:, 0],
                     _n * conc[:, 1],
                     meff,
                     meff,
                     n_e * self.mpi3tofmm3,
                     n_mu * self.mpi3tofmm3,
                     np.zeros(_n.shape),
                     np.zeros(_n.shape),
                     np.zeros(_n.shape),
                     self.mu_e * self.m_pi,
                     conc[:, 1]
        ]).transpose()
        tab = tabulate(table, ['mu_n [MeV]', 'P [MeV/fm^3]', 'E[MeV/fm^3]',
                               'n_b [1/fm^3]',
                               'n_n [1/fm^3]', 'n_p [1/fm^3]', 'm_n^*', 'm_p^*',
                               'n_e [1/fm^3]', 'n_mu [1/fm^3]', 'n_u [1/fm^3]',
                               'n_d [1/fm^3]', 'n_s [1/fm^3]',
                               'mu_e [MeV]', 'Y_p'], floatfmt='.6f', tablefmt='plain')
        f.write(tab)
        f.close()
        return table

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


    def mu(self, nrange=None):
        if nrange is None:
            nrange = self.nrange
        self.check()
        conc = self.rho
        mu = []
        for n in conc:
            mu.append([eos.mu(n, i+1, self.C) for i in range(self.n_baryon)])

        return arr(mu)

    def mu_gen(self, nrange_with_f):
        mu = []
        for n in nrange_with_f:
            mu.append([eos.mu(n, i+1, self.C) for i in range(self.n_baryon)])

        return arr(mu)

    def dumpMu(self):
        mu = self.m_pi*self.mu()
        tab = np.insert(mu, 0, self.nrange/self.n0, axis=1)
        np.savetxt(join(self.foldername, self.filenames['mu']), tab)

    def inspect_f(self):
        self.check()

        def swArr(n):
            out = eos.dArray(len(n))
            for i, _n in enumerate(n):
                out[i] = _n
            return out

        def func(f, n):
            arg = eos.dArray(1)
            out = eos.dArray(1)
            arg[0] = f
            params = eos.func_f_eq_params()

            params.n = n
            params.C = self.C
            params.dimN = self.n_baryon
            params.df = 1e-3
            eos.func_f_eq(arg, out, 1, 1, params)
            return out[0]
        frange = np.linspace(0, 1, 500)
        fig, ax = plt.subplots()
        l, = plt.plot(frange, list(map(
            lambda z: func(z, swArr(self.rho[0, 1:])), frange
        )))
        plt.subplots_adjust(left=0.25, bottom=0.25)
        axn = plt.axes([0.25, 0.1, 0.65, 0.03])
        sn = Slider(axn, 'N', 0, self.rho.shape[0]-1, valinit=0)

        def update(val):
            l.set_ydata(list(map(
            lambda z: func(z, swArr(self.rho[sn.val, 1:])), frange
        )))
            fig.canvas.draw_idle()
        sn.on_changed(update)
        plt.show()


class MassDriver():
    def __init__(self):
        self.N = None
        self.E = None
        self.P = None

    def setEos(self, N, E, P):
        self.N = N
        self.E = E
        self.P = P

class Nucleon(Wrapper):
    def __init__(self, C, basefolder_suffix=''):
        super(Nucleon, self).__init__(C, basefolder_suffix=basefolder_suffix)
        self.part_names=['n', 'p']
        self.C.Hyper = 0
        self.filenames['meff'] = 'meff.dat'
        self.filenames['mass_nocrust'] = 'masses_nocrust.dat'
        self.filenames['mass_crust'] = 'masses_crust.dat'
        self.filenames['eta'] = 'eta_ns.dat'
        self.filenames['vs'] = 'vs_ns.dat'
        self.filenames['Eparts'] = 'parts_ns.dat'
        self.filenames['Pparts'] = 'p_parts_ns.dat'
        self.filenames['uniparts'] = 'parts_uni_ns.dat'
        self.filenames['profiles'] = 'profiles_nucl.dat'
        self.filenames['grig'] = 'grig_nucl.dat'
        self.filenames['mu'] = 'mu_nucl.dat'


    def E(self, nrange, ret_f=False, f=0.):
        E, P, n = self.EPN(nrange=nrange)
        if ret_f:
            flist = self.rho[:, 0]
            return [E, flist]
        else:
            return [E]

    def dumpEos(self):
        pass

    def dumpMeff(self):
        self.check()
        tab = [[self.C.phi_n(i, self.C.Xs(i, z)*z)
                  for i in range(self.n_baryon)] for z in self.rho[:, 0]]
        data = np.insert(arr(tab), 0, self.nrange/self.n0, axis=1)
        table = tabulate(data, ['n/n0']+[self.part_names[i] for i in range(self.n_baryon)],
                         tablefmt='plain')
        with open(join(self.foldername, self.filenames['meff']), 'w') as f:
            f.write(table)



class Sym(Nucleon):
    def __init__(self, C, basefolder_suffix=''):
        super(Sym, self).__init__(C, basefolder_suffix=basefolder_suffix)
        self.C.Hyper = 0
        self.filenames.update(J='j_tilde.dat')
        self.filenames['eta'] = 'eta_sym.dat'
        self.filenames['vs'] = 'vs_sym.dat'
        self.filenames['Eparts'] = 'parts_sym.dat'
        self.filenames['Pparts'] = 'p_parts_sym.dat'
        self.filenames['uniparts'] = 'parts_uni_sym.dat'
        self.filenames['profiles'] = 'profiles_sym.dat'
        self.filenames['mu'] = 'mu_sym.dat'



    def reset(self):
        self._E, = self.E(self.nrange)
        self._P = self.P(self.nrange)
        self.set = 1

    def E(self, nrange, ret_f=False, f=0.):
        nlist = [[n / 2, n / 2] for n in nrange]
        dn = 1e-4
        dn = arr([dn, dn])
        return Wrapper.E_gen(self, nlist, ret_f, f, dn=dn)

    def P(self, nrange):
        nlist = []
        f = 0.
        for _n in nrange:
            f, = eos.f_eq(np.array([_n / 2, _n / 2]), np.array([f]), 1, self.C)
            nlist.append(np.array([f, _n / 2, _n / 2]))
        res = self.P_chem(nlist, leptons=False)
        return res

    def Jtilde(self, nrange=None):
        if nrange is None:
            nrange = self.nrange
        return arr([eos.J(z, self.C) for z in nrange])

    def dumpJ(self):
        tab = arr([self.nrange/self.n0, self.Jtilde(self.nrange)]).transpose()
        with open(join(self.foldername, self.filenames['J']), 'w') as f:
            f.write(tabulate(tab, ['n/n0', 'J\tilde[MeV]'], tablefmt='plain'))

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
                                         self.n0, dx=1e-1, n=2, order=9)




class Neutr(Nucleon):
    def __init__(self, C, basefolder_suffix=''):
        super(Neutr, self).__init__(C, basefolder_suffix=basefolder_suffix)
        self.C.Hyper = 0
        self.filenames['mass_crust'] = 'masses_crust_neutron.dat'
        self.filenames['mass_nocrust'] = 'masses_neutron.dat'
        self.filenames['meff'] = 'meff_ns.dat'
        self.filenames['eta'] = 'eta_n.dat'
        self.filenames['vs'] = 'vs_n.dat'
        self.filenames['Eparts'] = 'parts_neutr.dat'
        self.filenames['Pparts'] = 'p_parts_neutr.dat'
        self.filenames['uniparts'] = 'parts_uni_neutr.dat'
        self.filenames['profiles'] = 'profiles_neutr.dat'
        self.filenames['mu'] = 'mu_neutr.dat'


    def reset(self):
        self._E, = self.E(self.nrange)
        self._P = self.P(self.nrange)
        self.set = 1

    def E(self, nrange, ret_f=False, f=0.):
        nlist = [[n, 0] for n in nrange]
        dn = 1e-4
        dn = arr([dn, 0])
        return Wrapper.E_gen(self, nlist, ret_f, f, dn=dn)

    def P(self, nrange):
        nlist = []
        f = 0.
        for _n in nrange:
            f, = eos.f_eq(np.array([_n, 0]), np.array([f]), 1, self.C)
            nlist.append(np.array([f, _n, 0]))
        res = self.P_chem(nlist, leptons=False)
        return res




class Hyperon(Nucleon):
    def __init__(self, C, basefolder_suffix=''):
        super(Hyperon, self).__init__(C, basefolder_suffix=basefolder_suffix)
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
        self.filenames['grig'] = 'grig_hyper.dat'
        self.filenames['mu'] = 'mu_hyper.dat'
        self.part_names += ['L', 'Sm', 'S0', 'Sp', 'Xm', 'X0']

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
        n_e, n_mu = self.lepton_concentrations()
        self.dumpGrig(E, P, conc, n_e, n_mu)

        return tab

    def dumpEtap(self):
        frange = np.linspace(0, 1., 100)
        tab = arr([frange, [self.C.eta_p(z) for z in frange]]).transpose()
        with open(join(self.foldername, self.filenames['etap_f']), 'w') as f:
            f.write(tabulate(tab, ['f', 'eta_p(f)'], tablefmt='plain'))

    def dumpChi(self):
        pass


class HyperonPhi(Hyperon):
    def __init__(self, C, basefolder_suffix=''):
        super(HyperonPhi, self).__init__(C, basefolder_suffix=basefolder_suffix)
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
        self.filenames['grig'] = 'grig_hyper_phi.dat'
        self.filenames['mu'] = 'mu_hyper_phi.dat'






class HyperonPhiSigma(HyperonPhi):
    def __init__(self, C, basefolder_suffix=''):
        super(HyperonPhiSigma, self).__init__(C, basefolder_suffix=basefolder_suffix)
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
        self.filenames['grig'] = 'grig_hyper_phi_sigma.dat'
        self.filenames['mu'] = 'mu_hyper_phi_sigma.dat'

class DeltaBase(Wrapper):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.foldername = join(self.foldername, 'delta')
        if not os.path.exists(self.foldername):
            os.makedirs(self.foldername)
        if hasattr(self.C, 'setDeltaConstants'):
            self.C.setDeltaConstants(2, 0)
            self.part_names.append('D-')
            self.part_names.append('D0')
            self.part_names.append('D+')
            self.part_names.append('D++')
            self.n_baryon = 12

    def getSymm(self, n, lastx=0., lastf=0.):
        def eq(x, params):
            ### n = [n_p, n_D0]
            # x = [x['n_p'].value, x['n_D0'].value]
            n = params[0]
            n_ch = params[0]/2
            if x[0] < 0:
                return [100500., 0]
            n = arr([n/2 - x[0]/2, n/2 - x[0]/2, 0., 0., 0., 0., 0., 0., 0., x[0]/2, x[0]/2, 0.])
            f = eos.f_eq(n, arr([lastf]), 1, self.C)[0]
            n_in = np.insert(n, 0, f)
            # return np.sum(np.arr([
            #     eos.mu(1, n_in, self.C) - eos.mu(10, n_in, self.C),
            #     eos.mu(2, n_in, self.C) - eos.mu(11, n_in, self.C)
            # ])**2)
            res = arr([
                eos.mu(n_in, 1, self.C) - eos.mu(n_in, 10, self.C),
            ])
            # print(x, res)
            return [res, f]

        res = leastsq(lambda z: eq(z, [n])[0], [lastx])[0]
        # res = optimize.minimize(lambda z: eq(z, [n])[0], [lastx], bounds=[[0,None]], method='SLSQP').x
        if res > n:
            res = [0.]
        return [res] + eq(res, [n])

    def dumpDeltaSym(self, suffix=''):
        lastn = 0.
        lastf = 0.
        res = []
        for i, n in enumerate(self.nrange):
            _res = self.getSymm(n, lastx=lastn, lastf=lastf)
            print('res = ', _res)
            if (i > 1):
                lastn = _res[0][0]
            print('lastn =', lastn)
            lastf = _res[-1]
            res.append(_res)
        res = np.array(res)
        n_d = res[:, 0]
        # n_d[0][0] = [0.0]
        frac = []
        for i, n in enumerate(n_d):
            nd = n[0]
            nn = (self.nrange[i] - nd)/2
            frac.append([nn, nn, 0., 0., 0., 0., 0., 0., 0., nd/2, nd/2])
        frac = arr(frac)

        E = self.E_gen(frac)
        Ebind = self.m_pi * (E/self.nrange - self.C.M[0])
        n_w_f = np.insert(frac, 0, res[:, -1], axis=1)
        mu = self.mu_gen(n_w_f)
        P = self.P_chem(n_w_f)
        np.savetxt(join(self.foldername, 'delta_sym'+suffix+'.dat'),
                   arr([self.nrange/self.n0, E, P, [q[0]/self.nrange[i] for i,q in enumerate(n_d)], Ebind]).transpose())
        np.savetxt(join(self.foldername, 'delta_sym_mu'+suffix+'.dat'),
                   np.insert(mu, 0, self.nrange/self.n0, axis=1).transpose())



class Delta(DeltaBase, Hyperon):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        # Delta.__init__(self, C, basefolder_suffix='')

class DeltaPhi(DeltaBase, HyperonPhi):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)


class DeltaPhiSigma(DeltaBase, HyperonPhiSigma):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)

class Model(Wrapper):
    def __init__(self, C, K0=None, f0=None, J0=None, suffix=None, basefolder_suffix=''):
        if any([K0, f0, J0]):
            params = []
            names = []
            #a dirty workaround
            for p in [['K0', K0], ['J0', J0], ['f0', f0]]:
                if p[1] is not None:
                    names.append(p[0])
                    params.append(p[1])

            args = dict(list(zip(names, params)))
            print(args)
            wr = Wrapper_old(C())
            wr.solve(**args)

        super(Model, self).__init__(C, basefolder_suffix=basefolder_suffix)
        self.nucl = Nucleon(C, basefolder_suffix=basefolder_suffix)
        self.sym = Sym(C, basefolder_suffix=basefolder_suffix)
        self.neutr = Neutr(C, basefolder_suffix=basefolder_suffix)
        self.hyper = Hyperon(C, basefolder_suffix=basefolder_suffix)
        self.hyper_phi = HyperonPhi(C, basefolder_suffix=basefolder_suffix)
        self.hyper_phi_sigma = HyperonPhiSigma(C, basefolder_suffix=basefolder_suffix)
        self.delta = Delta(C, basefolder_suffix=basefolder_suffix)
        self.delta_phi = DeltaPhi(C, basefolder_suffix=basefolder_suffix)
        self.delta_phi_sigma = DeltaPhiSigma(C, basefolder_suffix=basefolder_suffix)
        self.children = [self.nucl, self.sym, self.neutr, self.hyper,
                         self.hyper_phi, self.hyper_phi_sigma,
                         self.delta, self.delta_phi]

        if any([K0, f0, J0]):
            C = wr.C
            for subwr in [self] + self.children:
                subwr.C.Cs = C.Cs
                subwr.C.Co = C.Co
                subwr.C.Cr = C.Cr
                subwr.C.b = C.b
                subwr.C.c = C.c
                subwr.C.f0 = C.f0

        # if suffix is not None:
        #     self.foldername = os.path.join(BASEFOLDER+basefolder_suffix,
        #                                    self.Ctype.__name__.strip('_'))
        #     if not os.path.exists(self.foldername):
        #         os.makedirs(self.foldername)
        #     for child in self.children:
        #         child.foldername = self.foldername
        self.filenames['eos'] = 'eos.dat'
        self.filenames.update(eta='eta_F.dat')
        self.filenames.update(props='props.dat')
        self.filenames.update(intro='intro.dat')


    def setDeltaConst(self, Xs, Xo, Xr, folder_suffix):
        for m in [self.delta, self.delta_phi, self.delta_phi_sigma]:
            m.C.setDeltaRho(Xr)
            m.C.setDeltaOmega(Xo)
            m.C.setDeltaSigma(Xs)
            m.foldername = self.foldername + folder_suffix
            if not os.path.exists(m.foldername):
                os.makedirs(m.foldername)

    def setParams(self, Cs, Co, Cr, b, c, f0):
        for m in [self]+self.children:
            m.C.Cs = Cs
            m.C.Co = Co
            m.C.Cr = Cr
            m.C.b = b
            m.C.c = c
            m.C.f0 = f0


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

        table = self.nucl.dumpGrig(E, P, conc, n_e, n_mu)
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
            s.dumpVs()
            s.dumpMassesCrust(npoints=10)
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
        print(pots)
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
            print(pots)
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

    def dumpPotentials(self):
        i = 0
        slist = []
        vlist = []
        f = [0.]
        for n in self.nrange:
            n_in = np.array([n/2, n/2])
            f = eos.f_eq(n_in, np.array(f), 1, self.C)
            pots = eos.potentials(np.insert(n_in, 0, f), 5, self.C)
            V = self.m_pi*self.C.X_o[i] * pots[2]
            S = self.m_pi*(-self.C.M[i]) * (1 - self.C.phi_n(i, self.C.X_s[i]*self.C.M[0]/self.C.M[i]*f[0]))
            slist.append(S)
            vlist.append(V)
        tab = np.array([self.nrange/self.n0, slist, vlist]).transpose()
        table = tabulate(tab, ['n/n_0', 'S [MeV]', 'V [MeV]'], tablefmt='plain')
        print(table)
        print(join(self.foldername, 'pot_sym.dat'))
        with open(join(self.foldername, 'pot_sym.dat'), 'w') as f:
            f.write(table)

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
        print(pots)
        V = -self.C.X_o[i] * pots[2]
        S =  (-self.C.M[i]) * (1 - self.C.phi_n(i, self.C.X_s[i]*self.C.M[0]/self.C.M[i]*f[0]))
        Uopt = lambda z: z + self.C.M[0] - sqrt((z + self.C.M[0] - V)**2 - S*(2 * self.C.M[i] + S))
        erange = np.linspace(-50./self.m_pi, 10., 100)
#         plt.plot(erange, self.m_pi*Uopt(erange))
#         plt.show()
        return self.m_pi*erange, self.m_pi*Uopt(erange)

import eosWrap as eos
# from lmfit.minimizer import minimize
# from lmfit.parameter import Parameters
from matplotlib.widgets import Slider, Button
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
from scipy.optimize import root, bisect
from multiprocessing import Queue, Process
from Wrapper import Wrapper as Wrapper_old
import six
from scipy.optimize import minimize
# import ipdb
import inspect
from copy import copy
workfolder = '/home/const/MEGA/'
import pdb
BASEFOLDER = join(workfolder, 'data2/')

def initBasefolder(fname):
    global BASEFOLDER
    BASEFOLDER = join(workfolder, fname) + '/'


class Wrapper(object):
    def __init__(self, Ctype, basefolder_suffix=''):
        self.Ctype = Ctype
        self.C = Ctype()
        self.foldername = os.path.join(BASEFOLDER + basefolder_suffix,
                 Ctype.__name__.strip('_'))
        if not os.path.exists(self.foldername):
            os.makedirs(self.foldername)
        self.m_pi = 135.
        self.m_e = 0.5 / self.m_pi
        self.m_mu = 105.0 / self.m_pi
        self.m_o = 783./self.m_pi
        self.n0 = self.C.n0
        self.set = False
        self.driver_set = False
        self.folder = None
        self.n_baryon = 2
        self.verbose = False
        self.mpi4_2_mevfm3 = self.m_pi * (self.m_pi / 197.33) ** 3
        self.mpi3tofmm3 = (self.m_pi / 197.33)**3
        self.nmin = 0.
        self.nmax = 8 * self.n0
        self.npoints = 800

        gamma = 1
        self.nrange = np.linspace(self.nmin, self.nmax**gamma, self.npoints,
                                  endpoint=False)**(1/gamma)
# example of setting filenames for a descendant of Wrapper
        self.filenames = {'mass_nocrust': None, 'eos': None,
                          'mass_crust': None, 'meff': None,
                          'eta': None, 'vs': None,
                          'Eparts': None,
                          'Pparts': None,
                          'uniparts': None,
                          'profiles': None,
                          'grig': None,
                          'mu': None,
                          'pf': None}

        self.part_names = ['n', 'p', 'Lambda', 'Sigma-',
                           'Sigma0', 'Sigma+', 'Xi-', 'Xi0', 'e', 'mu']

        self.md = MassDriver()

    def stepE(self, n, last, f, lenLast, iter, C, que):
        rho = eos.stepE(n, last, f, len(last), iter, C)
        que.put(rho)

    def getDuCrit(self):
        self.check()
        nc = 8 * self.n0
        solve = 0
        for i, _n in enumerate(self.nrange):
            if eos.p_f(self.rho[i, 1], 0.5) < eos.p_f(self.rho[i, 2], 0.5) + eos.p_f(self.n_e[i], 0.5):
                nc = _n
                solve = 1
                break
        # out = self.dumpMassesCrust(write=0)
        try:
            out = np.loadtxt(join(self.foldername, self.filenames['mass_crust']), skiprows=1)
        except FileNotFoundError:
            out = np.loadtxt(join(self.foldername, self.filenames['mass_crust'])+'_linear', skiprows=1)
        iM = interp1d(out[:, 0], out[:, 1])

        iPfn = interp1d(self.nrange, 
               [eos.p_f(r, 0.5) for r in self.rho[:, 1]])

        iPfp = interp1d(self.nrange, 
               [eos.p_f(r, 0.5) for r in self.rho[:, 2]])

        iPfe = interp1d(self.nrange, 
               [eos.p_f(r, 0.5) for r in self.n_e])
        
        from scipy.optimize import root
        if solve:
            root = root(lambda z: iPfn(z) - iPfp(z) - iPfe(z), nc, )
            nc = root.x
        Mc = iM(nc/self.n0)
        return [nc, Mc]



    def getParams(self, inv=0):
        try:
            if not inv:
                self.loadEos()
            else:
                self.loadEosInv()
                self.switch_inv()

        except FileNotFoundError:
            print('EoS file not found; Aborting.')
            return
        try:
            mass_data = np.loadtxt(join(self.foldername,
                                        self.filenames['mass_crust']+'_linear'), skiprows=1)
        except FileNotFoundError:
            print('Mass data file not found; Aborting.')
            return
        conc = self.concentrations()
        arg_crit = [np.where(c > 1e-6)[0][0] if np.where(c > 1e-6)[0].size > 0
        else len(c)-1 for c in conc.transpose()]
        n_crit = self.nrange[arg_crit]/self.n0

        n_cr = self.nrange[np.where(self.nc/self.nrange > 1e-6)[0][0]]/self.n0

        n_DU, m_DU = self.getDuCrit()
        iM = interp1d(mass_data[:, 0], mass_data[:, 1], kind='cubic')
        iR = interp1d(mass_data[:, 0], mass_data[:, 2], kind='cubic')
        
        nmax = minimize(lambda z: -iM(z), 7.).x[0]
        mmax = iM(nmax)
        rmax = iR(nmax)

        n15 = root(lambda z: iM(z)-1.5, [3.]).x[0]
        m15 = iM(n15)
        r15 = iR(n15)

        

        return n_crit, n_cr, n_DU/self.n0, m_DU, [nmax, mmax, rmax], [n15, m15, r15]
        

    def processMaxw(self, mu_init=None, show=0, branch_3=0, shift=0, verbose=0, mu=None):
        if mu is None:
            mu = np.nan_to_num(self.mu(branch_3=branch_3)[:, 0])
        if not branch_3:
            P = self._P
            E = self._E
            N = self.nrange
        else:
            P = self._P2
            E = self._E2
            N = self.nrange

        regions = []
        i_dP = np.where(np.diff(mu) < 0)[0]
        i0 = i_dP[0]
        region_lasti = 0
        for i, _n in enumerate(i_dP):
            if (_n - i0) > 1:
                regions.append(i_dP[region_lasti:i])
                region_lasti = i
                i0 = _n
            else:
                i0 = _n
        
        regions.append(i_dP[region_lasti:])
        
        i_break =  shift + np.where(np.diff(mu[shift:]) < 0)[0][-1]
        i_break_mu = shift + np.where(np.diff(mu[shift:]) < 0)[0][0]
        if verbose: print(i_break, i_break_mu)
        b1 = np.array([mu[:i_break_mu], P[:i_break_mu], E[:i_break_mu], N[:i_break_mu]]).transpose()
        b2 = np.array([mu[i_break+1:], P[i_break+1:], E[i_break+1:], N[i_break+1:]]).transpose()
        # for a in [mu[i_break+1:], P[i_break+1:], E[i_break+1:], N[i_break+1:]]:
            # print(a.shape, a.dtype)
        # print([mu[:i_break_mu], P[:i_break_mu], E[:i_break_mu], N[:i_break_mu]])
        # print([mu[i_break+1:], P[i_break+1:], E[i_break+1:], N[i_break+1:]])
        # return b1, b2
        # kind = 'cubic'
        kind = 'linear'
        # print(b1)
        # print(b2)
        # print(b1.shape, b2.shape)
        ip1 = interp1d(b1[:, 0], b1[:, 1], bounds_error=0, fill_value=0., kind=kind)
        ip2 = interp1d(b2[:, 0], b2[:, 1], bounds_error=0, fill_value=0., kind=kind)
        ie1 = interp1d(b1[:, 0], b1[:, 2], bounds_error=0, fill_value=0., kind=kind)
        ie2 = interp1d(b2[:, 0], b2[:, 2], bounds_error=0, fill_value=0., kind=kind)
        in1 = interp1d(b1[:, 0], b1[:, 3], bounds_error=0, fill_value=0., kind=kind)
        in2 = interp1d(b2[:, 0], b2[:, 3], bounds_error=0, fill_value=0., kind=kind)

        # print(b1[:i_break, 0], b2[:, 0])

        # mu_inter = np.intersect1d(b1[:, 0], b2[:, 0])
        mu_inter = np.linspace(min(b2[:, 0]), max(b1[:, 0]))
        # print(mu_inter)
        i_eq = np.argmin(abs(ip1(mu_inter) - ip2(mu_inter)))

        if verbose: print('i_eq = ', i_eq)
        if mu_init is None:
            # mu_init = [mu[i_break]]
            mu_init = [mu_inter[i_eq]]
            if verbose: print('mu_init = ', mu_init)
        if show:
            plt.plot(mu_inter, ip1(mu_inter))
            plt.plot(mu_inter, ip2(mu_inter))
            plt.plot(b2[:, 0], ip2(b2[:, 0]))
            plt.plot(b1[:, 0], ip1(b1[:, 0]))
            plt.show()
        # mu_eq = bisect(lambda z: ip1(z) - ip2(z), mu_init-1., mu_init+1.)
            plt.plot(b1[:, 0], b1[:, 1])
            plt.plot(b2[:, 0], b2[:, 1])
            plt.show()
        mu_eq = root(lambda z: ip1(z) - ip2(z), mu_init, tol=1e-6).x
        if verbose: print('mu_eq = ', mu_eq)
        P_eq = ip1(mu_eq)
        P_eq2 = ip2(mu_eq)
        if verbose: 
            print('P_eq = ', P_eq, 'P_eq2 = ', P_eq2)
            print('E1 = ', ie1(mu_eq))
            print('E2 = ', ie2(mu_eq))
        n1 = in1(mu_eq)
        n2 = in2(mu_eq)
        if verbose: 
            print('n1 = %.6f, n2 = %.6f' %(n1/self.n0, n2/self.n0))
            print('E1 = ', (ie1(mu_eq)/n1 - self.C.M[0]) * self.m_pi)
            print('E2 = ', (ie2(mu_eq)/n2 - self.C.M[0]) * self.m_pi)

        p1 = np.array([p for p in b1[:, 1] if p < P_eq])
        p2 = np.array([p for p in b2[:, 1] if p > P_eq])

        mu1 = np.array([mu for mu in b1[:, 0] if mu < mu_eq])
        mu2 = np.array([mu for mu in b2[:, 0] if mu > mu_eq])

        n1 = b1[:len(p1), 3]
        n2 = b2[(len(b2[:, 1]) - len(p2)):, 3]

        e1 = b1[:len(p1), 2]
        e2 = b2[(len(b2[:, 1]) - len(p2)):, 2]

        p_total = np.concatenate((p1, p2))
        n_total = np.concatenate((n1, n2))
        e_total = np.concatenate((e1, e2))
        mu_total = np.concatenate((mu1, mu2))

        self.nrange_maxw = n_total
        self._P_maxw = p_total
        self._E_maxw = e_total
        self.mu_maxw = mu_total



    def loadEos(self):
        raise NotImplementedError("Method should be overloaded")

    def reset(self, iterations=30, timeout=None):
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
                    self.set = 1
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
            rho[i] = np.insert(rho[i], 0, f)  # and again sp omitted
            mu_e.append(eos.mu(rho[i], 1, self.C) -
                        eos.mu(rho[i], 2, self.C))

        self.rho = np.ascontiguousarray(arr(rho))
        self.mu_e = np.ascontiguousarray(arr(mu_e))
        eparts = []
        _E = []
        for z in self.rho:
            epart = np.zeros((9), dtype='float')
            _E.append(eos.E(z, self.C))
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

    def dumpDensities(self):
        n_V = [sum([self.C.X_o[i] * r[i+1]/self.n0 for i in range(self.n_baryon)])
                for r in self.rho]
        n_I = [sum([self.C.X_r[i] * self.C.T[i] * r[i+1]/self.n0 for i in range(self.n_baryon)])
                for r in self.rho]
        n_S = [sum([self.C.X_p[i] * r[i+1]/self.n0 for i in range(self.n_baryon)])
                for r in self.rho]

        n_Q = [sum([self.C.Q[i] * r[i+1]/self.n0 for i in range(self.n_baryon)])
                for r in self.rho]

        np.savetxt(join(self.foldername, self.filenames['density']),
            np.array([self.nrange/self.n0, n_V, n_I, n_S, n_Q]).transpose())


    def getFnBranches(self):
        self.check()
        branches = [np.zeros(self.nrange.shape[0])]
        f_prev = 0.

        back_div = 2.
        num_space = 40
        num_roots = 1
        fmax = 1.
        for i, n in enumerate(self.rho[:, 1:]):
            fspace = np.linspace(f_prev/back_div, fmax, num_space)
            roots = []
            for f in fspace:
                root = eos.f_eq(arr(n), arr([f]), 1, self.C)[0]
                roots.append(root)
            print(n, roots)
            roots = sorted(set(roots))
            print(roots)
            _num_roots = len(set(roots))
            if _num_roots > num_roots:
                delta = _num_roots - num_roots
                num_roots = _num_roots
                for r in range(delta):
                    branches.append(np.zeros(self.nrange.shape[0]))
            for j, root in enumerate(roots):
                branches[j][i] = root

        return branches


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


    def dumpMassesCrust(self, nmin=0.4, nmax=None, npoints=100, write=True, fname=None, ret_frac=False):
        inter = 'linear'
        if nmax == None:
            nmax = 10.5*self.n0
        if self.needsMaxw() and not (hasattr(self, 'nrange_maxw')):
            raise ValueError('Maxwell construction inside dumpMassesCrust is deprecated, prepare the EoS explicitly instead!')
            # self.processMaxw()
            # if any(np.diff(self._P_maxw) < 0) or any(np.diff(self.nrange_maxw) < 0):
            #     raise ValueError('EoS still non-monotonous after Maxwell construction, check the solution!')
            # else:
            #     self.switch_maxw()
        
        out = self.stars_crust(nmin=nmin, nmax=nmax, npoints=npoints, ret_frac=ret_frac, inter=inter,
            show=0)
        out[0] /= self.n0
        if write:
            tab = arr(out).transpose()
            names = ['n/n0', 'M[M_sun]', 'R[km]', 'Mg_rect[M_sun]', 'Mg_trap[M_sun]']
            if ret_frac:
                names += self.part_names[1:self.n_baryon]
            table = tabulate(tab, names, tablefmt='plain')
            if fname is None:
                fname = self.filenames['mass_crust']+'_'+inter
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

    def stars_crust(self, ncut_crust=0.45, ncut_eos=0.7, inter='linear',
            n_stars=None, nmin=.6, nmax=5.0, npoints=50,
                    crust="crust.dat", show=False, crustName=None,
                    ret_frac=False, fasthyp=False, neutron=0, ret_i=0, 
                    force_reset=0, internal_kind=0):
        neos = self.npoints
        if nmax > max(self.nrange):
            nmax = max(self.nrange)
        if n_stars is not None:
            nmax = n_stars[-1]  # TODO too dirty workaround

        # print('in stars_crust')
        if not np.all([self.md.E, self.md.N, self.md.P]) or force_reset:
            self.check(nrange=self.nrange)
            E, N, P, e, finalE, finalN, finalP, n, p = self.setCrust(inter,
            ncut_crust, ncut_eos)
            self.md.setEos(finalN, finalE, finalP)
            
        else:
            finalE = self.md.E
            finalN = self.md.N
            finalP = self.md.P

        # return finalN, finalE, finalP 

        mpi2km = 5.7202e-5

        # finalN_low, finalE, finalP = self.md.rawN, self.md.rawE, self.md.rawP
        finalN_low, finalE, finalP = self.md.N, self.md.E, self.md.P
        # print(finalE)
        # pdb.set_trace()
        # return finalE, finalP, finalN_low*self.n0
        self.dr = eos.KVDriver(finalE, finalP, finalN_low*self.n0, internal_kind)

        if show:
            drE = arr([self.dr.EofN(z) for z in finalN*self.n0])
            drP = arr([self.dr.PofN(z) for z in finalN*self.n0])
            drDe = arr([
                derivative(self.dr.EofN, z, dx=1e-6, order=3) for z in finalN*self.n0
            ])
            drDe2 = arr([self.dr.dEofN(z) for z in finalN*self.n0])
            # np.savetxt('MKVOR_4EE.dat', np.array([finalN*self.n0, drE, drP]).transpose(), fmt='%.16e')
            print(drP)
            # lines = plt.plot(finalN, drE, n, arr(e)/self.m_pi**4, N, arr(E)/self.m_pi**4)
            lines = plt.plot(drE, drP, arr(e)/self.m_pi**4, arr(p)/self.m_pi**4, arr(E)/self.m_pi**4, arr(P)/self.m_pi**4)
            plt.show()

            lines = plt.plot(finalN, drE/finalN, n, arr(e)/n/self.m_pi**4, N, arr(E)/N/self.m_pi**4)
            # lines = plt.plot(finalN, finalP, finalN)
            plt.xlabel(r'$n/n_0$', fontsize=18)
            plt.xlim([0, 2])
            plt.ylabel(r'$P \, [MeV^4] $', fontsize=18)
            plt.legend(lines, ['interpolated', 'crust', 'RMF'], loc=0)
            plt.show()

            # np.savetxt('driver_crustKVOR.dat', np.array([finalN, drE, drP]).transpose(), fmt='%.8f')
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
            # with open(crustName, 'w') as f:
            #     f.write(table)

        if n_stars is None:
            n_stars = np.linspace(nmin, nmax, npoints)

        #interpolating particle fractions (except n)

        conc = self.concentrations()
        if ret_frac:
            inter_hyp = [interpolate.interp1d(self.nrange, conc[:, i])
                         for i in range(1, self.n_baryon)]

        #         plt.plot(self.n, map(lambda z: [f(z) for f in inter_hyp], self.n))
        #         plt.show()
        MR = []
        Mgrav = []
        str_fracs = []


        for _n in n_stars:
            mr = eos.star_crust2(_n, 3, self.dr, 1e-10)
            # print(_n/self.n0, mr)
            MR.append(mr)
            lastN = self.dr.getLastN(self.dr.nSize)[:-1]
            lastR = self.dr.getLastR(self.dr.nSize)[:-1]
            lastM = self.dr.getLastM(self.dr.nSize)[:-1]
            # try:
            dx = lastR[1] - lastR[0]
            # except IndexError:
            #     dx=1.

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

    def setCrust2(self, inter, ncut_crust, ncut_eos):
        E, P, N = self.EPN(self.nrange)
        P = P / self.mpi4_2_mevfm3
        N = N / self.n0
        # print(self.Ctype.__name__)
        np.savetxt('eos' + self.Ctype.__name__ + '.dat', arr([E, P, N]).transpose())
        E *= self.m_pi ** 4
        P *= self.m_pi ** 4
        e = [0.]
        p = [0.]
        n = [0.]
        # e = []
        # p = []
        # n = []
        with open("/home/const/workspace/swigEos/crust.dat", 'r') as f:
            for line in f:
                # print line
                _e, _p, _n = line.split()
                if float(_n) < ncut_crust:
                    e.append(float(_e))
                    p.append(float(_p))
                    n.append(float(_n))
        crust = np.loadtxt("/home/const/workspace/swigEos/crust.dat")
        crust[:, 0] /= self.m_pi ** 4
        crust[:, 1] /= self.m_pi ** 4
        np.savetxt('crust_export.dat', crust)
        n_eos = 5
        i_n_eos = np.argmin(abs(N - [ncut_eos for i in N]))
        plist = np.append(p[:], P[i_n_eos:(i_n_eos + n_eos)])
        elist = np.append(e[:], E[i_n_eos:(i_n_eos + n_eos)])
        nlist = np.append(n[:], N[i_n_eos:(i_n_eos + n_eos)])

        nraw = np.append(n[:], N[i_n_eos:])
        praw = np.append(p[:], P[i_n_eos:])
        eraw = np.append(e[:], E[i_n_eos:])
        conc = self.concentrations()[i_n_eos:]
        cr_conc = arr([[1.] + [0. for i in range(conc.shape[1]-1)] for n in e])
        # print(conc.shape, cr_conc.shape)
        res_conc = np.append(cr_conc, conc, axis=0).transpose()
        self.md.setRawEos(nraw, eraw/self.m_pi**4, praw/self.m_pi**4, res_conc)
        # print(i_n_eos)
        # print(P[i_n_eos:(i_n_eos + n_eos)])
        # print(N[i_n_eos:(i_n_eos + n_eos)])
        # exit()
        iP = interpolate.interp1d(nlist, plist, kind=inter)
        iE = interpolate.interp1d(nlist, elist, kind=inter)
        gamma = 1. / 4.
        iN = np.linspace(0, ncut_eos ** gamma, 1000)
        iN = iN ** (1. / gamma)
        crust_p = np.array(list(map(iP, iN)))
        crust_e = np.array(list(map(iE, iN)))
        np.savetxt("crust_dense.dat", np.array([crust_e/self.m_pi**4,
                                                crust_p/self.m_pi**4, iN]).transpose())
        #         finalE = np.append(crust_e, E[i_n_eos+n_eos:])
        #         finalP = np.append(crust_p, P[i_n_eos + n_eos:])
        #         finalN = np.append(iN, N[i_n_eos+n_eos:])
        finalE = np.append(crust_e, E[i_n_eos + n_eos:]) / self.m_pi ** 4
        finalP = np.append(crust_p, P[i_n_eos + n_eos:]) / self.m_pi ** 4
        finalN = np.append(iN, N[i_n_eos + n_eos:])
        finalE[0] = 0
        finalP[0] = 0
        return E, N, P, e, finalE, finalN, finalP, n, p

    def dumpPotentials(self):
        if not self.set:
            raise
        slist = []
        vlist = []
        rlist = []
        for r in self.rho:
            pots = eos.potentials(r, 5, self.C)
            V = [self.m_pi*self.C.X_o[i] * pots[2] for i in range(self.n_baryon)]
            S = [self.m_pi*(-self.C.M[i]) * (1 - self.C.phi_n(i, self.C.X_s[i]*self.C.M[0]/self.C.M[i]*r[0])) for i in range(self.n_baryon)]
            R = [self.m_pi*self.C.X_r[i]*self.C.T[i] * pots[3] for i in range(self.n_baryon)]
            slist.append(S)
            vlist.append(V)
            rlist.append(R)
        slist = np.array(slist)
        vlist = np.array(vlist)
        rlist = np.array(rlist)

        tS = np.array([self.nrange/self.n0] + slist.transpose().tolist()).transpose()
        tV = np.array([self.nrange/self.n0] + vlist.transpose().tolist()).transpose()
        tR = np.array([self.nrange/self.n0] + rlist.transpose().tolist()).transpose()

        tabS = tabulate(tS, ['n/n_0'] + ['S [' + pname + ']' for pname in self.part_names],
            tablefmt='plain')
        tabV = tabulate(tV, ['n/n_0'] + ['V [' + pname + ']' for pname in self.part_names],
            tablefmt='plain')
        tabI = tabulate(tR, ['n/n_0'] + ['I [' + pname + ']' for pname in self.part_names],
            tablefmt='plain',   )
        # return tab
        # print(table)
        # print(join(self.foldername, 'pot_sym.dat'))
        with open(join(self.foldername, self.filenames['S']), 'w') as f:
            f.write(tabS)
        with open(join(self.foldername, self.filenames['V']), 'w') as f:
            f.write(tabV)
        with open(join(self.foldername, self.filenames['I']), 'w') as f:
            f.write(tabI)
        return [arr(slist), arr(vlist)]


    def setCrust(self, inter, ncut_crust, ncut_eos):
        # E, P, N = self.EPN(self.nrange)
        if not self.needsMaxw():
            print('not switching to maxw')
            E, P, N = self.EPN(self.nrange)
        else:
            if not hasattr(self, 'nrange_maxw'):
                raise ValueError('Need to set up maxw first!')
            print('switching to maxw')
            E = np.copy(self._E_maxw)
            P = np.copy(self._P_maxw)
            N = np.copy(self.nrange_maxw)

        P = P / self.mpi4_2_mevfm3
        N = N / self.n0
        # print(self.Ctype.__name__)
        # np.savetxt('eos' + self.Ctype.__name__ + '.dat', arr([E, P, N]).transpose())
        E *= self.m_pi ** 4
        P *= self.m_pi ** 4
        e = [0.]
        p = [0.]
        n = [0.]
        # e = []
        # p = []
        # n = []
        with open("/home/const/Numerics/swigEos/crust.dat", 'r') as f:
            for line in f:
                # print line
                _e, _p, _n = line.split()
                if float(_n) < ncut_crust:
                    e.append(float(_e))
                    p.append(float(_p))
                    n.append(float(_n))
        crust = np.loadtxt("/home/const/Numerics/swigEos/crust.dat")
        crust[:, 0] /= self.m_pi ** 4
        crust[:, 1] /= self.m_pi ** 4
        # np.savetxt('crust_export.dat', crust)
        n_eos = 5
        i_n_eos = np.argmin(abs(N - [ncut_eos for i in N]))
        plist = np.append(p[:], P[i_n_eos:(i_n_eos + n_eos)])
        elist = np.append(e[:], E[i_n_eos:(i_n_eos + n_eos)])
        nlist = np.append(n[:], N[i_n_eos:(i_n_eos + n_eos)])

        nraw = np.append(n[:], N[i_n_eos:])
        praw = np.append(p[:], P[i_n_eos:])
        eraw = np.append(e[:], E[i_n_eos:])

        # iEP = interp1d(eraw/self.m_pi**4, praw/self.m_pi**4, kind='cubic')
        # erange = np.linspace(0., max(eraw/self.m_pi**4), 1000)
        # np.savetxt('pe_' + self.__repr__() + '.dat', arr([erange, iEP(erange)]).transpose())
        # np.savetxt('np_' + self.__repr__() + '.dat', arr([]))

        conc = self.concentrations()[i_n_eos:]
        cr_conc = arr([[1.] + [0. for i in range(conc.shape[1]-1)] for n in e])
        # print(conc.shape, cr_conc.shape)
        res_conc = np.append(cr_conc, conc, axis=0).transpose()
        self.md.setRawEos(nraw, eraw/self.m_pi**4, praw/self.m_pi**4, res_conc)
        # print(i_n_eos)
        # print(P[i_n_eos:(i_n_eos + n_eos)])
        # print(N[i_n_eos:(i_n_eos + n_eos)])
        # exit()
        iP = interpolate.interp1d(nlist, plist, kind=inter)
        iE = interpolate.interp1d(nlist, elist, kind=inter)
        gamma = 1. / 4.
        iN = np.linspace(0, ncut_eos ** gamma, 1000)
        iN = iN ** (1. / gamma)
        crust_p = np.nan_to_num(list(map(iP, iN)))
        crust_e = np.nan_to_num(list(map(iE, iN)))
        # print(crust_e, crust_p)
        # exit()
        # np.savetxt(self.__repr__() + "crust_dense.dat", np.array([crust_e/self.m_pi**4,
        #                                         crust_p/self.m_pi**4, iN]).transpose())
        #         finalE = np.append(crust_e, E[i_n_eos+n_eos:])
        #         finalP = np.append(crust_p, P[i_n_eos + n_eos:])
        #         finalN = np.append(iN, N[i_n_eos+n_eos:])
        finalE = np.append(crust_e, E[i_n_eos + n_eos:]) / self.m_pi ** 4
        finalP = np.append(crust_p, P[i_n_eos + n_eos:]) / self.m_pi ** 4
        finalN = np.append(iN, N[i_n_eos + n_eos:])
        finalE[0] = 0
        finalP[0] = 0
        # print(finalE)
        # exit()
        return E, N, P, e, finalE, finalN, finalP, n, p

    def dumpScalingsN(self):
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


    def E_gen(self, nlist, solve_f=1, ret_f=False, f=0., dn=None):
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
        print(ret_f)
        for n in nlist:
            epart = np.zeros(9, dtype='float')
            if not solve_f:
                ret_f = 0
                elist.append(eos._E(n, self.C))
            else:
                f = eos.f_eq(arr(n), arr([f]), 1, self.C)[0]
                flist.append(f)
                elist.append(eos._E(np.insert(n, 0, f), self.C))
            eparts.append(epart)
            self.Eparts = arr(eparts)
            if dn is not None:

                # dEparts.append(np.gradient(part, sum(nlist[1]) - sum(nlist[0])))
                de1 = np.zeros(9, dtype='float')
                de2 = np.zeros(9, dtype='float')
                eos._E(np.insert(n+dn/2, 0, f), self.C)
                eos._E(np.insert(n-dn/2, 0, f), self.C)
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
        if ret_f:
            E = out[0]
        else:
            E = out
        res = np.nan_to_num(self.m_pi * (E / nrange - self.C.M[0]))
        # if ret_f:
        #     return res, out[1]
        # else:
        #     return res
        if not ret_f:
            return res
        else:
            return [res, out[1]]

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
        iP = interp1d(n_eos, P)
        with open(join(self.foldername,    '%.3f'%(n_c/self.n0)+self.filenames['profiles']), 'w') as f:
            f.write(tabulate(arr([n/self.n0, r] + part + [mg] + [gm] + [iGamma(n), iP(n)]).transpose(), ['n/n0', 'r[km]']
                           + self.part_names[1:] + ['M', 'Mb', 'Gamma', 'P'], tablefmt='plain'))
            # f.write(tabulate(arr([n/self.n0, r] + [mg] + [gm] + [iGamma(n)]).transpose(), ['n/n0', 'r[km]', 'M[M_sun]']
            #                + ['gamma', 'Gamma'], tablefmt='plain', floatfmt='.3f'))


    def dumpGrig(self, E, P, conc, n_e, n_mu):
        f = open(join(self.foldername, self.filenames['grig']), 'w')
        _n = self.nrange * self.mpi3tofmm3
        mu_n = [self.m_pi * eos.mu(z, 1, self.C) for z in self.rho]
        meff = [self.C.phi_n(0, z) for z in self.rho[:, 0]]
        mH = [ [self.C.phi_n(0, self.C.X_s[i] * self.C.M[0]/self.C.M[i] * z) for z in self.rho[:, 0]]
               for i in range (2, self.n_baryon)]
        table = arr([mu_n, P, E * self.mpi4_2_mevfm3,
                     _n,
                     _n * conc[:, 0],
                     _n * conc[:, 1],
                     meff,
                     meff,
                     n_e * self.mpi3tofmm3,
                     n_mu * self.mpi3tofmm3]+
                    3 * [np.zeros(self.nrange.shape)]
                    +[
                     self.mu_e * self.m_pi,
                     conc[:, 1]
        ] + [_n*conc[:,i] for i in range(2, self.n_baryon)] + mH).transpose()
        tab = tabulate(table, ['mu_n [MeV]', 'P [MeV/fm^3]', 'E[MeV/fm^3]',
                               'n_b [1/fm^3]',
                               'n_n [1/fm^3]', 'n_p [1/fm^3]', 'm_n^*', 'm_p^*',
                               'n_e [1/fm^3]', 'n_mu [1/fm^3]', 'n_u [1/fm^3', 'n_d [1/fm^3',
                               'n_s [1/fm^3]'] +  [
                               'mu_e [MeV]', 'Y_p'] + ['Y_' + i for i in self.part_names[2:]]
                       + ['m_' + i for i in self.part_names[2:]], floatfmt='.6f', tablefmt='plain')
        f.write(tab)
        f.close()
        return table

    def dumpFortin(self):
        self.reset()
        E, P, n = self.EPN()
        conc = self.concentrations()
        table = np.array([n/self.n0 * 0.16, E * self.mpi4_2_mevfm3, P] +
                         [conc[:, i] for i in range(conc.shape[1])]).transpose()
        print(table)
        tab = tabulate(table, ['n[fm^-3]', 'E[MeV/fm^3]' ,'P[MeV/fm^3]'] + self.part_names, tablefmt='plain')
        with open(join(self.foldername,self.filenames['fortin']), 'w') as f:
            f.write(tab)

        # Output with crust
        self.setCrust( ncut_crust=0.45, ncut_eos=0.7, inter='cubic')

        table = arr([self.md.rawN * 0.16, self.md.rawE*self.mpi4_2_mevfm3,
                     self.md.rawP*self.mpi4_2_mevfm3] +
                     [self.md.rawConc[i] for i in range(self.n_baryon)]).transpose()
        print(self.md.rawConc)
        print(table)
        tab = tabulate(table, ['n[fm^-3]', 'E[MeV/fm^3]' ,'P[MeV/fm^3]'] + self.part_names,
                       tablefmt='plain')
        with open(join(self.foldername,'crust_'+self.filenames['fortin']), 'w') as f:
            f.write(tab)
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
        ne_list = arr(ne_list)/self.nrange
        nmu_list = arr(nmu_list)/self.nrange
        du_list = arr(du_list)
        if ret_du:
            return [ne_list, nmu_list, du_list]
        else:
            return np.array([ne_list, nmu_list]).transpose()
            # return [ne_list, nmu_list]

    def mu(self, nrange=None, branch_3=0):
        if nrange is None:
            nrange = self.nrange

        if branch_3:
            conc = self.rho2

        else:
            conc = self.rho
            self.check()
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
        np.savetxt(join(self.foldername, self.filenames['mu']), tab, fmt='%.6f')

    def needsMaxw(self):
        self.check()
        return any(np.diff(self._P)<0)

    def dumpPf(self, write=1):
        self.check()
        pf = arr([[eos.p_f(n, 2*self.C.S[i] + 1) for i, n in enumerate(r[1:])] for r in self.rho[0:]])
        ne, nmu = self.lepton_concentrations().transpose()
        pfe = arr([eos.p_f(n, 2.) for n in ne])
        pfmu = arr([eos.p_f(n, 2.) for n in nmu])

        pf = np.insert(pf, pf.shape[1], pfe, axis=1)
        pf = np.insert(pf, pf.shape[1], pfmu, axis=1)

        tab = np.insert(self.m_pi*pf, 0, self.nrange/self.n0, axis=1)
        if write:
            np.savetxt(join(self.foldername, self.filenames['pf']), tab, fmt='%.6f')
        return self.m_pi*pf

    def get_feq(self, n):
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
        IofN = lambda z: int(interp1d(self.nrange, range(len(self.nrange)))(z))
        print(IofN(n))
        func = list(map(
            lambda z: func(z, swArr(self.rho[IofN(n), 1:])), frange
        ))

        return frange, func

    def func_f(self, f, n, mu):
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

    def swArr(self, n):
        out = eos.dArray(len(n))
        for i, _n in enumerate(n):
            out[i] = _n
        return out

    def inspect_f(self, rho=None, mu_e=None, show=True, i=None):
        self.check()
        if rho is None:
            rho = self.rho
        if mu_e is None:
            mu_e = self.mu_e
        if i is None:
            show = 1

        if show:
            frange = np.linspace(0, 1, 500)
            fig, ax = plt.subplots()

            # func = list(map(
            #     lambda z: func(z, swArr(self.rho[0, 1:])), frange
            # ))
            if i is None:
                i = 1
            l, = ax.plot(frange, list(map(
                lambda z: self.func_f(z, self.swArr(rho[i, 1:]), mu_e[i]), frange
            )))

            ax.plot(frange, [0. for f in frange])

            plt.subplots_adjust(left=0.25, bottom=0.25)
            axn = plt.axes([0.25, 0.1, 0.65, 0.03])
            ax_button = plt.axes([0.25, 0.05, 0.65, 0.03])
            btn = Button(ax_button, "Get values")
            sn = Slider(axn, 'N', 0, rho.shape[0]-1, valinit=i)

            def update(val):
                l.set_ydata(list(map(
                lambda z: self.func_f(z, self.swArr(rho[sn.val, 1:]), mu_e[sn.val]), frange
            )))
                fig.canvas.draw_idle()

            def get_values(event):
                global res
                print("get values!")
                values = l.get_ydata()
                x = l.get_xdata()
                res = [x, values]
                # print(res)
            res = []
            btn.on_clicked(get_values)
            sn.on_changed(update)
            ax.set_ylim([-20, 7])
            plt.show()
        else:
            frange = np.linspace(0, 1, 500)

            res = arr(
                [self.func_f(z, self.swArr(rho[i, 1:]), mu_e[i]) for z in frange]
                )
            return frange, res


        # res = [None, None]

        # print(res)
        return res



    def inspect_f_indep(self):
        self.check()

        frange = np.linspace(0, 1, 500)
        fig, ax = plt.subplots()

        # func = list(map(
        #     lambda z: func(z, swArr(self.rho[0, 1:])), frange
        # ))

        l, = ax.plot(frange, list(map(
            lambda z: self.func_f(z, self.swArr(self.rho[0, 1:]), self.mu_e[0]), frange
        )))

        ax.plot(frange, [0. for f in frange])

        plt.subplots_adjust(left=0.25, bottom=0.25)
        axn = plt.axes([0.25, 0.1, 0.65, 0.03])
        sn = Slider(axn, 'N', 0, self.rho.shape[0]-1, valinit=0)

        axi = plt.axes([0.25, 0.1, 0.65, 0.03])
        si = Slider(axi, 'i', 0, self.rho.shape[0]-1, valinit=0)

        def update_i(val):
            l.set_ydata(list(map(
            lambda z: self.func_f(z, self.swArr(self.rho[sn.val, 1:]), self.mu_e[sn.val]), frange
        )))

        def update_n(val):
            pass
            fig.canvas.draw_idle()
        si.on_changed(update_i)
        ax.set_ylim([-20, 7])
        plt.show()

    def checkMultiSol(self, nrange=None, rho=None, mu=None, npoints=10, ret_3br=0, break_lost=0):
        if nrange is None:
            nrange = self.nrange
        if rho is None:
            rho = self.rho
        # if mu is None:
        #     mu = self.mu

        n_check = nrange[::len(nrange)/npoints]
        i_check = range(0, len(nrange), len(nrange)//npoints)
        eqs = []
        for i in i_check:
            frange, res = self.inspect_f(show=0, rho=rho, mu_e=mu, i=i)
            eqs.append(res)

        # n_out = [] #densities at which more that 1 solution occurs
        # f_out = [] #approximate values of the roots for these densities
        out = []
        flag = 0
        sch_prev = 0
        for i_n, res in enumerate(eqs):
            sign_changes = 0
            r_prev = res[0]
            f_res = []
            for i, r in enumerate(res):
                if r*r_prev < 0:
                    sign_changes += 1
                    f_res.append(frange[i])
                    r_prev = r
            if sign_changes > 1 - flag:
                out.append([sign_changes, i_check[i_n]] + f_res)
                if ret_3br:
                    flag = 1
            if sign_changes < sch_prev:
                if break_lost:
                    return out
            sch_prev = sign_changes
        return out




class MassDriver():
    def __init__(self):
        self.N = None
        self.E = None
        self.P = None
        self.rawN = None
        self.rawE = None
        self.rawP = None
        self.rawConc = None

    def setEos(self, N, E, P):
        self.N = N
        self.E = E
        self.P = P

    def setRawEos(self, N, E, P, conc):
        self.rawN = N
        self.rawP = P
        self.rawE = E
        self.rawConc = conc

class Nucleon(Wrapper):
    def __init__(self, C, basefolder_suffix=''):
        self.name = 'nucl'
        super(Nucleon, self).__init__(C, basefolder_suffix=basefolder_suffix)
        self.part_names=['n', 'p']
        self.C.Hyper = 0
        self.filenames['meff'] = 'meff.dat'
        self.filenames['mass_nocrust'] = 'masses_nocrust.dat'
        self.filenames['mass_crust'] = 'masses_crust.dat'
        self.filenames['eta'] = 'eta_nucl.dat'
        self.filenames['vs'] = 'vs_ns.dat'
        self.filenames['Eparts'] = 'parts_ns.dat'
        self.filenames['Pparts'] = 'p_parts_ns.dat'
        self.filenames['uniparts'] = 'parts_uni_ns.dat'
        self.filenames['profiles'] = 'profiles_nucl.dat'
        self.filenames['grig'] = 'grig_nucl.dat'
        self.filenames['mu'] = 'mu_nucl.dat'
        self.filenames['fortin'] = 'fortin_nucl.dat'
        self.filenames['pf'] = 'pf_nucl.dat'
        self.filenames['eos'] = 'nucl.dat'
        self.filenames['S'] = 'S_nucl.dat'
        self.filenames['V'] = 'V_nucl.dat'
        self.filenames['I'] = 'I_nucl.dat'
        self.filenames['density'] = 'dens_nucl.dat'





    def loadEos(self):
        try:
            eos_tab = np.loadtxt(join(self.foldername, self.filenames['eos']), skiprows=1)
        except:
            raise
        self.nrange = eos_tab[:, 0] * self.n0
        self._E = eos_tab[:, 1]
        self._P = eos_tab[:, 2]
        self.set = 1
        rhos = self.nrange*eos_tab[:, 3 : 3 + self.n_baryon].transpose()
        # print(rhos)
        self.rho = np.insert(rhos.transpose(), 0, eos_tab[:,3+self.n_baryon], axis=1)

        mu_e = []
        n_e = []
        n_mu = []
        for r in self.rho:
            # print(r)
            # print(eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C))
            _mue = eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C)
            mu_e.append(_mue)
            if _mue > self.m_e:
                n_e.append((_mue**2 - self.m_e**2)**(3./2) / (3 * pi **2))
            else:
                n_e.append(0.)
            if _mue > self.m_mu:
                n_mu.append((_mue**2 - self.m_mu**2)**(3./2) / (3*pi**2))
            else:
                n_mu.append(0.)

        self.n_e = np.array(n_e)
        self.n_mu = np.array(n_mu)

        # print(self.n_e +self.n_mu - self.rho[:,2])
        # exit()

        self.mu_e = np.array(mu_e)
        self.mu_c = np.array(mu_e)

        self.xDUp = ((3*pi**2 * self.rho[:, 1])**(1./3) - (3*pi**2 * self.n_e)**(1./3))**(3.) / (3*pi**2 * self.nrange)

    def E(self, nrange, ret_f=False, f=0.):
        E, P, n = self.EPN(nrange=nrange)
        if ret_f:
            flist = self.rho[:, 0]
            return [E, flist]
        else:
            return [E]

    def dumpEos(self):
        E, P, n = self.EPN()
        tab = [self.nrange/self.n0, E, P]
        conc = self.concentrations()
        tab += [conc[:,i] for i in range(self.n_baryon)]
        tab += [self.rho[:,0]]
        tab += [self.m_pi*(E/n - self.C.M[0])]
        table = tabulate(arr(tab).transpose(), ['n/n0', 'E_NS[m_\pi^4]',
                                                'P_NS[m_\pi^4]'] +
        self.part_names[:self.n_baryon] + ['f', 'Ebind'], tablefmt='plain'
        )
        with open(join(self.foldername, self.filenames['eos']), 'w') as f:
            f.write(table)
        n_e, n_mu = self.lepton_concentrations().transpose()
        self.dumpGrig(E, P, conc, n_e, n_mu)


    def dumpMeff(self):
        print([self.C.Xs(i, 2*self.n0) for i in range(8)])
        self.check()
        tab = [[self.C.phi_n(i, self.C.Xs(i, z)*self.C.M[0]/self.C.M[i]*z)
                  for i in range(self.n_baryon)] for z in self.rho[:, 0]]
        data = np.insert(arr(tab), 0, self.nrange/self.n0, axis=1)
        table = tabulate(data, ['n/n0']+[self.part_names[i] for i in range(self.n_baryon)],
                         tablefmt='plain')
        with open(join(self.foldername, self.filenames['meff']), 'w') as f:
            f.write(table)


    def get_upper_branch(self, alpha=0.75):
        self.loadEos()
        out = self.checkMultiSol()
        if out == []:
            print("No other branches found")
            return
        f_init = np.array([out[-1][-1]])
        i_start = len(self.nrange) - 1
        # Third solution with selfconsistent concentrations usually appears
        # earlier than it was with "old" concentrations => factor alpha
        if alpha > 1:
            alpha = 1
        i_end = int(alpha * out[0][1])
        print(i_end)
        n_tosolve = self.nrange[i_start:i_end:-1]
        n_init = self.rho[i_start]
        init = n_init[2:]

        # Solving and populating concentrations
        flist = []
        rlist = []
        mu_e = []
        # self.nc2 = []
        for i_solve, _n in enumerate(n_tosolve):
        #     print(_n/m.n0)
            res = eos.stepE(_n, init, f_init, len(init), 30, self.C)
            n_S = np.array([_n - res[0], res[0]])
            # self.nc2.append(res[-1])
            f = eos.f_eq(n_S, f_init, 1, self.C)
        #     f = eos.f_eq(n_S, f_init, 1, m.C)
        #     print(res, f)
            rlist.append([_n - res[0], res[0]])
            flist.append(f[0])
            f_init = f

        rlist = np.array(rlist)
        flist = np.array(flist)
        # self.nc2 = np.array(self.nc2)
        self.rho_upper = np.insert(rlist, 0, flist.transpose(), axis=1)
        mu_e =np.array([eos.mu(r, 0, self.C) - eos.mu(r, 1, self.C) for r in self.rho_upper])
        self.mu_e_upper = mu_e
        self.n_upper = n_tosolve

        # print(self..checkMultiSol)

        # Stripping the new concentrations for the actual solution

    def refine_upper(self):
        sols = self.checkMultiSol(nrange=self.n_upper, rho=self.rho_upper,
            mu=self.mu_e_upper, npoints=len(self.n_upper), ret_3br=0)
        print(sols)
        i_break = sols[-1][1]
        print(i_break)
        self.rho_upper2 = self.rho_upper[:i_break]
        self.mu_e_upper2 = self.mu_e_upper[:i_break]
        self.n_upper2 = self.n_upper[:i_break]
        # self.nc2 = self.nc2[:i_break]
        # self.E_upper = np.array([eos.E_rho(n, self.mu_e_upper2[i], self.C) for i, n in enumerate(self.rho_upper2)])
        self._E2 = self._E.copy()
        self.rho2 = self.rho.copy()
        self.mu_e2 = self.mu_e.copy()
        self.E_upper = self.E_gen(self.rho_upper2, solve_f=0)
        self._P2 = self._P.copy()
        # self.nclist2 = self.nc.copy()
        self.P_upper = self.P_chem(self.rho_upper2)
        shift = len(self._E2) - len(self.E_upper)
        flag = 0 #VERY DIRTY!!!!!!!!!!!!!!!!!!!! assumes that third solution
                                                #doesn't disappear and is energetically
                                                #favorable for all densities!
        for i, e in enumerate(self.E_upper[::-1]):
            if self._E2[shift + i] > e:
                flag = 1
            if flag:
                self._E2[shift + i] = e
                self.rho2[shift + i] = self.rho_upper2[::-1][i]
                self.mu_e2[shift + i] = self.mu_e_upper2[::-1][i]
                self._P2[shift + i] = self.P_upper[::-1][i]
                # self.nclist2[shift + i] = self.nc2[::-1][i]


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
        self.filenames['eos'] = 'sym.dat'



    def reset(self):
        self._E, flist = self.E(self.nrange, ret_f=1)
        self._P = self.P(self.nrange)
        nlist = arr([ [n/2, n/2] for n in self.nrange]).transpose()
        # f = arr([0.])
        # flist = []
        # for n in nlist:
        #     f = eos.f_eq(n, f, 1, self.C)
        #     flist.append(f)
        # flist = arr(flist)
        print(nlist.shape, flist.shape)
        self.rho = np.insert(nlist, 0, flist, axis=0).transpose()
        self.mu_e = arr([0. for n in self.nrange])
        print(self.rho.shape)
        self.set = 1

    def E(self, nrange, ret_f=False, f=0.):
        nlist = [[n / 2, n / 2] for n in nrange]
        dn = 1e-4
        dn = arr([dn, dn])
        return Wrapper.E_gen(self, nlist, solve_f=1, ret_f=ret_f, f=f, dn=dn)

    def P(self, nrange):
        nlist = []
        f = 0.
        for _n in nrange:
            f, = eos.f_eq(np.array([_n / 2, _n / 2]), np.array([f]), 1, self.C)
            nlist.append(np.array([f, _n / 2, _n / 2]))
        res = self.P_chem(nlist, leptons=False)
        return res

    def Jtilde(self, nrange=None, f=None):
        if nrange is None:
            nrange = self.nrange
        return arr([eos.J(z, self.C) for z in nrange])

    def dumpJ(self):
        self.check()
        res = [[n/self.n0, eos.J(n, self.C, f)] for n, f in zip(self.nrange, self.rho[:, 0])]
        tab = arr(
                  res
                  )
        with open(join(self.foldername, self.filenames['J']), 'w') as f:
            f.write(tabulate(tab, ['n/n0', 'J\tilde[MeV]'], tablefmt='plain'))

        # return res

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
        self.filenames['eos'] = 'neutr.dat'
        self.filenames['mass_crust'] = 'masses_crust_neutron.dat'
        self.filenames['mass_nocrust'] = 'masses_neutron.dat'
        self.filenames['meff'] = 'meff_neutr.dat'
        self.filenames['eta'] = 'eta_neutr.dat'
        self.filenames['vs'] = 'vs_neutr.dat'
        self.filenames['Eparts'] = 'parts_neutr.dat'
        self.filenames['Pparts'] = 'p_parts_neutr.dat'
        self.filenames['uniparts'] = 'parts_uni_neutr.dat'
        self.filenames['profiles'] = 'profiles_neutr.dat'
        self.filenames['mu'] = 'mu_neutr.dat'


    def reset(self):
        self._E, flist = self.E(self.nrange, ret_f=1)
        self._P = self.P(self.nrange)
        self.rho = np.array([[f, n, 0.] for n, f in zip(self.nrange, flist)])
        self.mu_e = np.array([-1. for n in self.nrange])
        self.set = 1

    def E(self, nrange, ret_f=False, f=0.):
        nlist = [[n, 0] for n in nrange]
        dn = 1e-4
        dn = arr([dn, 0])
        return Wrapper.E_gen(self, nlist, solve_f=1, ret_f=ret_f, f=f, dn=dn)

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
        self.name = 'hyper'
        self.C.Hyper = 1
        self.n_baryon = 8
        self.C.phi_meson = 0
        self.C.hyper_phi_sigma = 1
        self.C.hyper_sigma_kind = 1
        self.filenames['mass_crust'] = 'mass_hyper.dat'
        self.filenames['mass_nocrust'] = None
        self.filenames['eos'] = 'hyper.dat'
        self.filenames.update(etap_f='etap_f.dat')
        self.filenames.update(etap_n='etap_n.dat')
        self.filenames['meff'] = 'meff_hyper.dat'
        self.filenames['vs'] = 'vs_hyper.dat'
        self.filenames['grig'] = 'grig_hyper.dat'
        self.filenames['mu'] = 'mu_hyper.dat'
        self.filenames['fortin'] = 'fortin_hyper.dat'
        self.filenames['pf'] = 'pf_hyper.dat'
        self.filenames['eta'] = 'eta_hyper.dat'
        self.filenames['S'] = 'S_hyper.dat'
        self.filenames['V'] = 'V_hyper.dat'
        self.filenames['I'] = 'I_hyper.dat'
        self.filenames['density'] = 'dens_hyper.dat'


        self.part_names += ['L', 'Sm', 'S0', 'Sp', 'Xm', 'X0']

    def dumpEos(self):
        E, P, n = self.EPN()
        tab = [self.nrange/self.n0, E, P]
        conc = self.concentrations()
        tab += [conc[:,i] for i in range(self.n_baryon)]
        tab += [self.rho[:,0]]
        tab += [self.m_pi*(E / self.nrange - self.C.M[0])]
        table = tabulate(arr(tab).transpose(), ['n/n0', 'E_NS[m_\pi^4]',
                                                'P_NS[m_\pi^4]'] +
        self.part_names[:self.n_baryon] + ['f'] + ['Ebind'], tablefmt='plain'
        )
        with open(join(self.foldername, self.filenames['eos']), 'w') as f:
            f.write(table)
        n_e, n_mu = self.lepton_concentrations().transpose()
        self.dumpGrig(E, P, conc, n_e, n_mu)

        # if hasattr(self, 'status'):
        #     np.savetxt(join(self.foldername, self.filenames['eos'])+'_status', self.status)

        # if hasattr(self, 'warnings'):
        #     np.save(join(self.foldername, self.filenames['eos'])+'_warnings', self.status)

        # if hasattr(self, 'warnings_inv'):
        #     np.save(join(self.foldername, self.filenames['eos'])+'_warnings_inv', self.status)


        return tab

    def dumpEtap(self):
        frange = np.linspace(0, 1., 100)
        tab = arr([frange, [self.C.eta_p(z) for z in frange]]).transpose()
        with open(join(self.foldername, self.filenames['etap_f']), 'w') as f:
            f.write(tabulate(tab, ['f', 'eta_p(f)'], tablefmt='plain'))

    def dumpChi(self):
        pass

    def loadEos(self):
        try:
            eos_tab = np.loadtxt(join(self.foldername, self.filenames['eos']), skiprows=1)
        except:
            return
        self.nrange = eos_tab[:, 0] * self.n0
        self._E = eos_tab[:, 1]
        self._P = eos_tab[:, 2]
        self.set = 1
        rhos = self.nrange*eos_tab[:, 3 : 3 + self.n_baryon].transpose()
        # print(rhos)
        self.rho = np.insert(rhos.transpose(), 0, eos_tab[:, 3 + self.n_baryon], axis=1)

        mu_e = []
        n_e = []
        n_mu = []
        for r in self.rho:
            # print(r)
            # print(eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C))
            _mue = eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C)
            mu_e.append(_mue)
            if _mue > self.m_e:
                n_e.append((_mue**2 - self.m_e**2)**(3./2) / (3 * pi **2))
            else:
                n_e.append(0.)
            if _mue > self.m_mu:
                n_mu.append((_mue**2 - self.m_mu**2)**(3./2) / (3*pi**2))
            else:
                n_mu.append(0.)

        self.n_e = np.array(n_e)
        self.n_mu = np.array(n_mu)

        # print(self.n_e +self.n_mu - self.rho[:,2])
        # exit()

        self.mu_e = np.array(mu_e)

        self.xDUp = ((3*pi**2 * self.rho[:, 1])**(1./3) - (3*pi**2 * self.n_e)**(1./3))**(3.) / (3*pi**2 * self.nrange)



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
        self.filenames['fortin'] = 'fortin_hyper_phi.dat'
        self.filenames['pf'] = 'pf_hyper_phi.dat'
        self.filenames['eta'] = 'eta_hyper_phi.dat'
        self.filenames['S'] = 'S_hyper_phi.dat'
        self.filenames['V'] = 'V_hyper_phi.dat'
        self.filenames['I'] = 'I_hyper_phi.dat'
        self.filenames['density'] = 'dens_hyper_phi.dat'


class HyperonPhi2(Hyperon):
    def __init__(self, C, basefolder_suffix=''):
        super(HyperonPhi2, self).__init__(C, basefolder_suffix=basefolder_suffix)
        self.C.Hyper = 1
        self.C.phi_meson = 1
        self.C.phi_kind = 0
        self.C.hyper_sigma_kind = 1
        self.filenames['mass_crust'] = 'mass_hyper_phi2.dat'
        self.filenames['mass_nocrust'] = None
        self.filenames['eos'] = 'hyper_phi2.dat'
        self.filenames['etap_f'] = 'etap_phi2_f.dat'
        self.filenames.update(etap_n='etap_phi2_n.dat')
        self.filenames['meff'] = 'meff_hyper_phi2.dat'
        self.filenames['vs'] = 'vs_hyper_phi2.dat'
        self.filenames['grig'] = 'grig_hyper_phi2.dat'
        self.filenames['mu'] = 'mu_hyper_phi2.dat'
        self.filenames['fortin'] = 'fortin_hyper_phi2.dat'
        self.filenames['pf'] = 'pf_hyper_phi2.dat'
        self.filenames['eta'] = 'eta_hyper_phi2.dat'
        self.filenames['potentials'] = 'pot_hyper_phi2.dat'





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
        self.filenames['fortin'] = 'fortin_hyper_phi_sigma.dat'
        self.filenames['pf'] = 'pf_hyper_phi_sigma.dat'
        self.filenames['eta'] = 'eta_hyper_phi_sigma.dat'
        self.filenames['S'] = 'S_hyper_phi_sigma.dat'
        self.filenames['V'] = 'V_hyper_phi_sigma.dat'
        self.filenames['I'] = 'I_hyper_phi_sigma.dat'
        self.filenames['density'] = 'dens_hyper_phi_sigma.dat'


class DeltaBase(Wrapper):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.foldername = join(self.foldername, 'Delta')
        if not hasattr(self, 'basefoldername'):
            self.basefoldername = self.foldername
        if not os.path.exists(self.foldername):
            os.makedirs(self.foldername)
        if hasattr(self.C, 'setDeltaConstants'):
            self.C.setDeltaConstants(2, 0)
            self.part_names.append('D-')
            self.part_names.append('D0')
            self.part_names.append('D+')
            self.part_names.append('D++')
            self.n_baryon = 12
        self.filenames['var_u'] = 'var_u.dat'


    def loadEos(self):
        try:
            eos_tab = np.loadtxt(join(self.foldername, self.filenames['eos']), skiprows=1)
            self.nrange = eos_tab[:, 0] * self.n0
            self._E = eos_tab[:, 1]
            self._P = eos_tab[:, 2]
            self.set = 1
            rhos = self.nrange*eos_tab[:, 3 : 3 + self.n_baryon].transpose() ##-2 due to Ebind in the file
            # print(rhos)
            self.rho = np.insert(rhos.transpose(), 0, eos_tab[:, 3 + self.n_baryon], axis=1)

            # mu_e = []
            n_e = []
            n_mu = []
            grig = np.loadtxt(join(self.foldername, self.filenames['grig']), skiprows=1)
            mu_e = grig[:, 13]/self.m_pi
            # print(mu_e)
            # exit()

            for i, r in enumerate(self.rho):
                # print(r)
                # print(eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C))
                # _mue = eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C)
                _mue = mu_e[i]
                if _mue > self.m_e:
                    n_e.append((_mue**2 - self.m_e**2)**(3./2) / (3 * pi **2))
                else:
                    n_e.append(0.)
                if _mue > self.m_mu:
                    n_mu.append((_mue**2 - self.m_mu**2)**(3./2) / (3*pi**2))
                else:
                    n_mu.append(0.)

            self.n_e = np.array(n_e)
            self.n_mu = np.array(n_mu)

            # print(self.n_e +self.n_mu - self.rho[:,2])
            # exit()

            self.mu_e = np.array(mu_e)

            self.xDUp = ((3*pi**2 * self.rho[:, 1])**(1./3) - (3*pi**2 * self.n_e)**(1./3))**(3.) / (3*pi**2 * self.nrange)
            self.set = 1

        except FileNotFoundError:
            raise
            # self.dumpEos()

    def getSymm(self, n, lastx=0., lastf=0.):
        def eq(x, n):
            if x[0] < 0:
                return [100500., 0]
            n = arr([n/2 - x[0]/2, n/2 - x[0]/2, 0., 0., 0., 0., 0., 0., x[0]/4, x[0]/4, x[0]/4, x[0]/4])
            f = eos.f_eq(n, arr([lastf]), 1, self.C)[0]
            n_in = np.insert(n, 0, f)
            res = [
                eos.mu(n_in, 1, self.C) - eos.mu(n_in, 10, self.C),
            ]

            return res + [f]
        res = leastsq(lambda z: eq(z, n)[0], [lastx], ftol=1e-16)[0].tolist()
        p_range = np.linspace(0, n, 100)
        if res[0] > 1e-5:
            pass
            # plt.plot(p_range, list(map(lambda z: eq([z], n)[0], p_range)))
            # plt.show()
        if res[0] > n:
            res = [0.]
        return res + eq(res, n)

    def dumpDeltaSym(self, suffix='', write=1):
        lastn = 0.
        lastf = 0.
        res = []
        for i, n in enumerate(self.nrange):
            _res = self.getSymm(n, lastx=lastn, lastf=lastf)
            print('res = ', _res)
            if (i > 1):
                lastn = _res[0]
            print('lastn =', lastn)
            lastf = _res[-1]
            res.append(_res)
        res = np.array(res)
        n_d = res[:, 0]
        frac = []
        for i, n in enumerate(n_d):
            nd = n
            nn = (self.nrange[i] - nd)/2
            frac.append([nn, nn, 0., 0., 0., 0., 0., 0., nd/4, nd/4, nd/4, nd/4])
        frac = arr(frac)

        E = self.E_gen(frac)
        self._E = np.copy(E)
        Eparts = self.Eparts
        Ebind = self.m_pi * (E/self.nrange - self.C.M[0])
        n_w_f = np.insert(frac, 0, res[:, -1], axis=1)
        self.rho = np.copy(n_w_f)
        mu = self.mu_gen(n_w_f)
        P = self.P_chem(n_w_f)
        self._P = np.copy(P)
        self.set = 1
        if write:
            np.savetxt(join(self.foldername, 'delta_sym'+suffix+'.dat'),
                       arr([self.nrange/self.n0, E, P, [q/self.nrange[i] for i,q in enumerate(n_d)], Ebind]).transpose())
            np.savetxt(join(self.foldername, 'delta_sym_mu'+suffix+'.dat'),
                       np.insert(self.m_pi*mu, 0, self.nrange/self.n0, axis=1))
            tab = [[self.C.phi_n(i, self.C.Xs(i, z)*self.C.M[0]/self.C.M[i]*z)
                for i in range(self.n_baryon)] for z in res[:, -1]]
            data = np.insert(arr(tab), 0, self.nrange/self.n0, axis=1)
            table = tabulate(data, ['n/n0']+[self.part_names[i] for i in range(self.n_baryon)],
                             tablefmt='plain')
            with open(join(self.foldername, 'meff_deltasym.dat'), 'w') as f:
                f.write(table)
        return frac, E, Ebind, Eparts, res[:,-1]


class DeltaSym(DeltaBase, Sym):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.foldername = join(self.foldername, 'DeltaSym')
        if not os.path.exists(self.foldername):
            os.makedirs(self.foldername)
        self.filenames['eos'] = 'eos_sym.dat'

    def eq(self, x, n, lastf=0.):
            if x[0] < 0:
                return [100500., 0]
            n = arr([n/2 - x[0]/2, n/2 - x[0]/2, 0., 0., 0., 0., 0., 0., x[0]/4, x[0]/4, x[0]/4, x[0]/4])
            f = eos.f_eq(n, arr([lastf]), 1, self.C)[0]
            n_in = np.insert(n, 0, f)
            res = [
                eos.mu(n_in, 1, self.C) - eos.mu(n_in, 10, self.C),
            ]

            return res + [f]

    def eq2(self, x, n, lastf=0.):
        pass

    def checkDeltaMultiSol(self, nrange=None, npoints=20, break_lost=0):
        # if mu is None:
        #     mu = self.mu
        if nrange is None:
            nrange = self.nrange
        nd_res = []
        out = []
        i_check = range(0, len(nrange), len(nrange)//npoints)
        for i in i_check:
            nd, res, flist = self.get_fn(nrange[i], f_init=0.5)
            # n_out = [] #densities at which more that 1 solution occurs
            # f_out = [] #approximate values of the roots for these densities
            sch_prev = 0
            sign_changes = 0
            r_prev = res[0]
            for k, r in enumerate(res):
                if r*r_prev < 0:
                    sign_changes += 1
                    print(i, sign_changes, nd[k])
                    nd_res.append(nd[k])
                    r_prev = r

            if sign_changes > 1:
                out.append([sign_changes, i, nrange[i]] + nd_res)
                return out

            if sign_changes < sch_prev:
                if break_lost:
                    return out
            sch_prev = sign_changes
        return out


    def get_fn(self, n, f_init=0.):
        fprev = f_init
        # global fprev
        # print(n)
        res = []
        nrange = np.linspace(0, n, 500)
        flist = []
        for _n in nrange:
            out, f = self.eq([_n], n, lastf=fprev)
            fprev = f
            flist.append(f)
            res.append(out)
        flist = arr(flist)
        return [nrange, np.array(res), flist]

    def inspect_eq(self):
        self.check()


        fig, ax = plt.subplots()

        # func = list(map(
        #     lambda z: func(z, swArr(self.rho[0, 1:])), frange
        # ))

        n = max(self.nrange)
        res0 = self.get_fn(0.01)
        l, = ax.plot(res0[0]/0.01, res0[1])

        plt.subplots_adjust(left=0.25, bottom=0.25)
        axn = plt.axes([0.25, 0.1, 0.65, 0.03])
        sn = Slider(axn, 'N', 0, len(self.nrange), valinit=1.)

        def update(val):
            print('val=', val)
            _n = self.nrange[int(val)]
            res = self.get_fn(4*_n, f_init=self.rho[int(val), 0])
            l.set_ydata(res[1])
            l.set_xdata(res[0]/_n)
            fig.canvas.draw_idle()

        sn.on_changed(update)
        ax.set_ylim([-5, 5])
        plt.show()

    def reset(self, iterations=30, timeout=100, stopnc=0):
        print([self.C.X_s[i] for i in range(12)])

        lastx = 0.
        lastf = 0.
        deltas = []
        flist = []
        rho = []
        for n in self.nrange:
            res = leastsq(lambda z: self.eq(z, n, lastf)[0], [lastx], ftol=1e-16)[0].tolist()
            lastx = res[0]
            if lastx > n:
                lastx = 0.
            lastf = self.eq([lastx], n, lastf)[1]
            if stopnc:
                if lastx > 1e-5:
                    return n
            flist.append(lastf)
            deltas.append(lastx)
            rho.append([lastf, (n-lastx)/2, (n-lastx)/2, 0., 0., 0., 0., 0., 0.,
                        lastx/4, lastx/4, lastx/4, lastx/4])
            print(n, lastx, lastf)

        self.rho = arr(rho)
        self._E = self.E_gen(self.rho[:, 1:])
        self._P = self.P_chem(self.rho)
        self.set = 1

    def E(self, nrange=None, ret_f=0):
        if nrange is None:
            nrange = self.nrange
        self.check(nrange=nrange)
        return self._E

    def checkEq(self, n, ndinit, finit):
        nrange = np.linspace(0, n, 1000)
        eqlist = arr([self.eq([x], n, lastf=finit)[0] for x in nrange])
        plt.plot(nrange, eqlist)
        plt.plot(nrange, [0 for i in nrange])
        plt.show()

    def dumpEos(self, nmax=None, npoints=None, write=True):
        self.check()
        E, P, n = self.EPN()
        frange = self.rho[:, 0]
        conc = self.concentrations()
        table = [n/self.n0, frange, self.Ebind(self.nrange), P, conc[:, 0]*2, conc[:, 10]*4, E]
        tab = tabulate(arr(table).transpose(), ['n/n0', 'E', 'P', 'N', 'Delta'], tablefmt='plain')
        with open(join(self.foldername, self.filenames['eos']), 'w') as f:
            f.write(tab)

    def loadEos(self):
        data = np.loadtxt(join(self.foldername, self.filenames['eos']), skiprows=1)
        self.nrange = self.n0 * data[:, 0]
        self._E = self.nrange * (data[:, 2]/self.m_pi + self.C.M[0])
        self._P = data[:, 3]# / self.mpi4_2_mevfm3
        z = np.array([0. for i in data[:, 1]])
        self.rho =(self.nrange * np.array([np.nan_to_num(data[:, 1]/self.nrange), data[:, 4] / 2, data[:, 4] / 2,
                             z, z, z, z, z, z, data[:, 5] / 4, data[:, 5] / 4,
                             data[:, 5] / 4, data[:, 5] / 4])).transpose()
        self.mu_e = np.array([0. for n in self.nrange])
        # print(self.rho.shape)
        # print(data.shape)
        # exit()
        self.set = 1
        return




    def dumpNc(self, urange):
        res = []
        for U in urange:
            xs = self.getDeltaXs(U)
            self.C.setDeltaSigma(arr([xs for i in range(4)]))
            nc = self.reset(stopnc=1)
            if nc is None:
                break
            res.append(nc/self.n0)

        np.savetxt(join(self.foldername, 'nc.dat'), arr([xsrange[:len(res)], res]).transpose(), fmt='%.6f')

class DeltaAsym(DeltaSym):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.foldername = join(self.foldername, 'DeltaAsym')
        if not os.path.exists(self.foldername):
            os.makedirs(self.foldername)
        self.filenames['eos'] = 'eos_asym.dat'
        self.asym = 0.

    def eq(self, x, n, lastf=0.):
            if x[0] < 0:
                return [100500., 0]
            n_n = (0.5 + self.asym) * n - x[0]/2
            n_p = (0.5 - self.asym) * n - x[0]/2
            n = arr([n_n, n_p, 0., 0., 0., 0., 0., 0., x[0]/4, x[0]/4, x[0]/4, x[0]/4])
            f = eos.f_eq(n, arr([lastf]), 1, self.C)[0]
            n_in = np.insert(n, 0, f)
            res = [
                eos.mu(n_in, 1, self.C) - eos.mu(n_in, 10, self.C),
            ]

            return res + [f]

    def reset(self, iterations=30, timeout=100, stopnc=0):
        print([self.C.X_s[i] for i in range(12)])

        lastx = 0.
        lastf = 0.
        deltas = []
        flist = []
        rho = []
        for n in self.nrange:
            res = leastsq(lambda z: self.eq(z, n, lastf)[0], [lastx], ftol=1e-16)[0].tolist()
            lastx = res[0]
            if lastx > n:
                lastx = 0.
            lastf = self.eq([lastx], n, lastf)[1]
            if stopnc:
                if lastx > 1e-5:
                    return n
            flist.append(lastf)
            deltas.append(lastx)
            n_n = (0.5 + self.asym) * n - lastx/2
            n_p = (0.5 - self.asym) * n - lastx/2
            rho.append([lastf, n_n, n_p, 0., 0., 0., 0., 0., 0.,
                        lastx/4, lastx/4, lastx/4, lastx/4])
            print(n, lastx, lastf)

        self.rho = arr(rho)
        self._E = self.E_gen(self.rho[:, 1:])
        self._P = self.P_chem(self.rho)
        self.set = 1



class Delta(DeltaBase, Hyperon):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        # Delta.__init__(self, C, basefolder_suffix='')



class Rho(Wrapper):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.nc = np.array([])
        self.foldername = join(self.foldername, 'rcond')
        self.basefoldername = self.foldername
        if not os.path.exists(self.foldername):
            os.makedirs(self.foldername)
        self.filenames['rcond'] = 'rcond.dat'
        self.filenames['n_rho'] = 'n_rho.dat'
        self.C.chi_r_prime = 0

    def set_eos_names(self, suffix):
        #keylist = ['eos', 'grig', 'rcond', 'n_rho', 'density']
        if not hasattr(self, 'olddict'):
            self.olddict = copy(self.filenames)
        for k in self.filenames.keys():
            try:
                self.filenames[k] = self.olddict[k] + suffix
            except TypeError:
                pass

    def getEparts():
        Eparts = []
        for r, mu_e in zip(self.rho, self.mu_e):
            E, out = eos.E_rho_parts(r, mu_e, self.C, len(self.part_names) +  6)

            Eparts.append(out)

        return np.array(Eparts)


    def switch_inv(self):
        self.mu_e = self.mu_e_inv
        self.rho = self.rho_inv
        self._P = self._P_inv
        self._E = self._E_inv
        self.nrange = self.nrange_inv
        self.n_e = self.n_e_inv
        self.n_mu = self.n_mu_inv
        self.nc = self.nc_inv

    def switch_maxw(self):
        self.nrange = self.nrange_maxw
        self._P = self._P_maxw
        self._E = self._E_maxw
        if self.rho is None:
            self.rho = self.rho_inv

    def E_NF(self, n, f, init, iterations=30):
        res = eos.stepE_rho_nf(n, f, init, len(init) + 1,
                            iterations, init[-1], self.C)
        n_n = n - sum(res[:-2])
        nc = res[-1]
        mu_e = res[-2]
        n_in = np.concatenate(([f, n_n], res[:-2]))
        return res, self.Efull(rlist = [n_in], nclist=[nc], mu_list=[mu_e])

    def regen_mue(self):
        mu = self.mu()
        mu_e = mu[:, 0] - mu[:, 1]
        n_e = []
        n_mu = []
        for i, r in enumerate(self.rho):
            # print(r)
            # print(eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C))
            # _mue = eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C)
            # mu_e.append(_mue)
            _mue = mu_e[i]
            if _mue > self.m_e:
                n_e.append((_mue**2 - self.m_e**2)**(3./2) / (3 * pi **2))
            else:
                n_e.append(0.)
            if _mue > self.m_mu:
                n_mu.append((_mue**2 - self.m_mu**2)**(3./2) / (3*pi**2))
            else:
                n_mu.append(0.)

        self.n_e = np.array(n_e)
        self.n_mu = np.array(n_mu)
        self.mu_e = mu_e


    def dumpNR(self):
        self.check()
        nr = []
        mn = self.C.M[0]
        C = self.C
        for i, r in enumerate(self.rho):
            _f = r[0]
            mu_e = self.mu_e[i]
            nr.append((2 * mn**2 * C.m_rho / C.Cr) * C.eta_r(_f)**0.5 * C.phi_n(0, _f)**2 /
             np.float64(C.chi_prime(_f)) * (1 - mu_e / (C.m_rho * np.float64(C.phi_n(0, _f)))))

        nr = arr(nr)

        out = np.array([self.nrange/self.n0, nr/self.n0]).transpose()

        np.savetxt(join(self.foldername, self.filenames['n_rho']), out)


    def func_f(self, f, n, mu_c):
        arg = eos.dArray(1)
        out = eos.dArray(1)
        arg[0] = f
        params = eos.func_f_eq_params_rho()

        params.n = n
        params.C = self.C
        params.dimN = self.n_baryon
        params.df = 1e-3
        params.mu_c = mu_c

        eos.func_f_eq_rho(arg, out, 1, 1, params)
        return out[0]

    def mu(self, nrange=None, branch_3=0, inv=0):
        if nrange is None:
            nrange = self.nrange

        if not inv:
            self.check()

        if inv:
            conc = self.rho_inv
            mu_e = self.mu_e_inv
        else:
            if branch_3:
                conc = self.rho2
                mu_e = self.mu_e
            else:
                conc = self.rho
                mu_e = self.mu_e


        mu = []
        for j, n in enumerate(conc):
            mu.append([eos.mu_rho(n, i+1, mu_e[j], self.C) for i in range(self.n_baryon)])
        return arr(mu)

    def getFnBranches(self):
        self.check()
        branches = [np.zeros(self.nrange.shape[0])]
        f_prev = 0.

        back_div = 2.
        num_space = 10
        num_roots = 1
        fmax = 1.
        for i, n in enumerate(self.rho[:, 1:]):
            fspace = np.linspace(f_prev/back_div, fmax, num_space)
            roots = []
            for f in fspace:
                root = eos.f_eq_rho(arr(n), arr([f]), 1, self.mu_e[i], self.C)[0]
                roots.append(root)
            print(n, roots)
            roots = sorted(set(roots))
            print(roots)
            _num_roots = len(set(roots))
            if _num_roots > num_roots:
                delta = _num_roots - num_roots
                num_roots = _num_roots
                for r in range(delta):
                    branches.append(np.zeros(self.nrange.shape[0]))
            for j, root in enumerate(roots):
                branches[j][i] = root

        return branches

    def mu_deriv(self, nrange=None):
        if nrange is None:
            nrange = self.nrange
        self.check()
        conc = self.rho
        mu = []
        for j, n in enumerate(conc):
            mu.append([eos.mu_rho(n, i+1, self.mu_e[j], self.C) for i in range(self.n_baryon)])
        return arr(mu)

    def dumpEos(self):
        super().dumpEos()

        rcond = arr([self.nrange/self.n0, self.nc/self.nrange]).transpose()
        np.savetxt(join(self.foldername, self.filenames['rcond'])
                   ,rcond, fmt='%.8f')


    def loadEos(self):
        print('rho loadEos')
        super().loadEos()

        rcond = np.loadtxt(join(self.foldername, self.filenames['rcond']))
        self.nc = self.nrange * rcond[:, 1]

        mu_e = []
        n_e = []
        n_mu = []

        grig = np.loadtxt(join(self.foldername, self.filenames['grig']),
            skiprows=1)
        mu_e = grig[:, 13]/self.m_pi
        # print(mu_e)
        # exit()

        n_e = []
        n_mu = []
        for i, r in enumerate(self.rho):
            # print(r)
            # print(eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C))
            # _mue = eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C)
            # mu_e.append(_mue)
            _mue = mu_e[i]
            if _mue > self.m_e:
                n_e.append((_mue**2 - self.m_e**2)**(3./2) / (3 * pi **2))
            else:
                n_e.append(0.)
            if _mue > self.m_mu:
                n_mu.append((_mue**2 - self.m_mu**2)**(3./2) / (3*pi**2))
            else:
                n_mu.append(0.)

        self.n_e = np.array(n_e)
        self.n_mu = np.array(n_mu)

        # print(self.n_e +self.n_mu - self.rho[:,2])
        # exit()

        self.mu_e = np.array(mu_e)
        self.mu_c = np.array(mu_e)
        self.n_e = np.array(n_e)
        self.n_mu = np.array(n_mu)

        # print(self.n_e +self.n_mu - self.rho[:,2])
        # exit()

        self.mu_e = np.array(mu_e)
        self.mu_c = np.array(mu_e)


    def reset(self, iterations=30, timeout=None, kind=1):
        """Calculates the EoS stored in the wrapper. Has to be overriden in
        Symm, Neutron.

        Calculates the particle composition, scalar field, energy density,
        pressure, ... for the given parameter set self.C. Wrapper.set is
        set to True after s succesful calculation.

        Note: \sigma^* is not supported in the current version.

        """
        print("Resetting " + self.__repr__() + " for model " + self.Ctype.__name__.strip('_'))

        rho = []

        mu_e = []
        mu_c = []
        nc = []
        if kind == 1:
            stepper = eos.stepE_rho
            init = arr([0. for i in range(self.n_baryon - 1)] + [0.])

        f = arr([0.])  # sp omitted
        self.status = []
        self.warnings = []
        for i, _n in enumerate(self.nrange):
            if timeout is None:
                # print('init= ', init)
                out, info = eos.stepE_rho_withf(_n, init, f, len(init)+2, iterations, 0., self.C, 10+len(init))
                # print('init -> ', init)
            else:
                raise NotImplementedError

            fun = info[-len(init):]
            info = info[:-len(init)]
            info_dict = {'init':init, 'i': i, 'f': f, 'info': info, 'out': out, 'fun':fun}
            self.status.append(info_dict)
            
            if int(info[6]) != 6:
                self.warnings.append(info_dict)
     
                
            if i % (len(self.nrange) / 20) == 0:
                print('.', end=' ')

            rho.append(out.copy()[:-3]) ## 3 reserved for f, mu_c and n_c
            rho[i] = np.insert(rho[i], 0, _n - np.sum(rho[i]))
            f = np.array([out[-3]])
            if self.verbose:
                pass  # TODO: some verbose output
            # rho contains all scalar fields as well
            rho[i] = np.insert(rho[i], 0, f)  #and again sp omitted
            mu_e.append(out[-2])
            mu_c.append(out[-2])
            nc.append(out.copy()[-1])
            init = np.delete(out[:-1], -2)

        self.mu_c = np.array(mu_c)
        self.nc = np.nan_to_num(np.ascontiguousarray(arr(nc)))
        self.rho = np.ascontiguousarray(arr(rho))
        self.mu_e = np.ascontiguousarray(arr(mu_e))
        eparts = []
        # _E = []
        # for i, z in enumerate(self.rho):
        #     eitem, epart = self.Efull(z, self.mu_e[i])
        #     _E.append(eitem)
        #     eparts.append(epart)
        # self._E = np.array(_E)
        self._E = self.Efull()
        self.Eparts = arr(eparts)
        dEparts = []
        for part in self.Eparts.transpose():
            dEparts.append(np.gradient(part, self.nrange[1]-self.nrange[0]))
        dEparts = arr(dEparts)
        # p1 = self.nrange * dEparts
        # self.Pparts = p1.transpose() - self.Eparts
        # self._E = np.array(map(lambda z: eos.E(z, self.C), self.rho))
        self._E = np.ascontiguousarray(self._E)
        self._P = np.ascontiguousarray(self.P_chem())
        self.set = True


    def reset_from_init(self, n_init, f_init, init_vals, npoints=100, nmax=8., iterations=30, timeout=None, kind=1):
        print("Resetting from init " + self.__repr__() + " for model " + self.Ctype.__name__.strip('_'))
        rho = []
        mu_e = []
        mu_c = []
        nc = []
        self.status = []
        self.warnings = []
        nmax *= self.n0

        stepper = eos.stepE_rho
        init = arr([0. for i in range(self.n_baryon - 1)] + [0.])

        n_tosolve = n_init
        if nmax > n_init[-1]:
            n_tosolve = np.append(n_init[:-1], np.linspace(n_init[-1], nmax,
                                                      npoints))
            init_vals = np.append(init_vals[:-1], np.array([init_vals[-1]]*npoints),
                                  axis=0)

            f_init = np.append(f_init[:-1], np.array([f_init[-1]]*npoints), axis=0)

            print(n_tosolve.shape)
            print(init_vals.shape)



        for i, res in enumerate(zip(n_tosolve, f_init, init_vals)):
            _n, f, init = res
            f = np.array([f])
            if timeout is None:
                # print('before ', init)
                out, info = eos.stepE_rho_withf(_n, init, f, len(init)+2, iterations, 0., self.C, 10 + len(init))
                # print('after ', init)
            else:
                raise NotImplementedError()
                # print init
            fun = info[-len(init):]
            info = info[:-len(init)]
            info_dict = {'init':init, 'i': i, 'f': f, 'info': info, 'out': out, 'fun':fun}
            self.status.append(info_dict)
            
            if int(info[6]) != 6:
                self.warnings.append(info_dict)
            
            # print(init)
            # if i % (len(self.nrange) / 20) == 0:
            #     print('.', end=' ')


            rho.append(out.copy()[:-3]) ## 3 reserved for f, mu_c and n_c
            rho[i] = np.insert(rho[i], 0, _n - np.sum(rho[i]))
            f = np.array([out[-3]])
            if self.verbose:
                pass  # TODO: some verbose output
            # rho contains all scalar fields as well
            rho[i] = np.insert(rho[i], 0, f)  #and again sp omitted
            mu_e.append(out[-2])
            mu_c.append(out[-2])
            nc.append(out.copy()[-1])
            init = np.delete(out[:-1], -2)
            if i >= len(n_init) and i < len(n_tosolve)-1:
                init_vals[i+1] = init
                f_init[i+1] = f[0]

        self.nrange = n_tosolve
        self.mu_c = np.array(mu_c)
        self.nc = np.nan_to_num(np.ascontiguousarray(arr(nc)))
        self.rho = np.ascontiguousarray(arr(rho))
        self.mu_e = np.ascontiguousarray(arr(mu_e))

        n_e = []
        n_mu = []
        for i, r in enumerate(self.rho):
            # print(r)
            # print(eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C))
            # _mue = eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C)
            # mu_e.append(_mue)
            _mue = self.mu_e[i]
            if _mue > self.m_e:
                n_e.append((_mue**2 - self.m_e**2)**(3./2) / (3 * pi **2))
            else:
                n_e.append(0.)
            if _mue > self.m_mu:
                n_mu.append((_mue**2 - self.m_mu**2)**(3./2) / (3*pi**2))
            else:
                n_mu.append(0.)

        self.n_e = np.array(n_e)
        self.n_mu = np.array(n_mu)
        eparts = []
        # _E = []
        # for i, z in enumerate(self.rho):
        #     eitem, epart = self.Efull(z, self.mu_e[i])
        #     _E.append(eitem)
        #     eparts.append(epart)
        # self._E = np.array(_E)
        self._E = self.Efull()
        self.Eparts = arr(eparts)
        dEparts = []
        for part in self.Eparts.transpose():
            dEparts.append(np.gradient(part, self.nrange[1]-self.nrange[0]))
        dEparts = arr(dEparts)
        # p1 = self.nrange * dEparts
        # self.Pparts = p1.transpose() - self.Eparts
        # self._E = np.array(map(lambda z: eos.E(z, self.C), self.rho))
        self._E = np.ascontiguousarray(self._E)
        self._P = np.ascontiguousarray(self.P_chem())
        self.set = True

    def get_init_inv(self):
        if not hasattr(self, 'nrange_inv'):
            try:
                self.loadEosInv()
            except:
                raise NotImplementedError

        n_init = self.nrange_inv
        f_init = self.rho_inv[:, 0]
        init = np.insert(self.rho_inv[:, 2:], self.rho_inv.shape[1]-2, self.mu_e_inv[:], axis=1)

        return init, f_init, n_init

    def reset_from_inv(self, nmax_inv = 4., nmax=8., npoints_inv=1000, npoints_add=1000, iterations=100):
        self.reset_inv(nmax=nmax_inv, npoints=npoints_inv, iterations=iterations)
        n_init = self.nrange_inv
        f_init = self.rho_inv[:, 0]
        init = np.insert(self.rho_inv[:, 2:], self.rho_inv.shape[1]-2, self.mu_e_inv[:], axis=1)
        self.reset_from_init(n_init, f_init=f_init, init_vals=init, nmax=nmax, npoints=npoints_add, iterations=iterations)

    def reset_from_branches(self, branches, npoints=1000, nmax=4., iterations=30, kind='linear'):
        nmin = min([min(b['n']) for b in branches])

        # remove branches with only 1 point
        print([len(b['n']) for b in branches])
        for i, b in enumerate(branches):
            if len(b['n']) < 2:
                branches.pop(i)
        print([len(b['n']) for b in branches])
        _nmax = max([max(b['n']) for b in branches])
        if _nmax < nmax:
            nmax = _nmax

        #density range for which the `exact` solution will be found
        nrange = np.linspace(nmin, nmax, npoints)

        # set up the interpolators
        i_b_f = [interp1d(b['n'], b['f']) for b in branches]
        i_b_E = [interp1d(b['n'], b['E']) for b in branches]
        i_b_mu = [interp1d(b['n'], b['mu_e']) for b in branches]
        i_b_conc = [
            [interp1d(b['n'], conc) for conc in b['conc'].transpose()]
            for b in branches
        ]

        init = []
        f_init = []
        for n in nrange:
            Es = []
            for i, b in enumerate(branches):
                try:
                    E = i_b_E[i](n)
                except ValueError:
                    E = np.inf
                Es.append(E)
            i_min = np.argmin(Es)
            # print(n, Es, i_min)
            f_init.append(i_b_f[i_min](n))
            init.append(np.array([i_conc(n) for i_conc in i_b_conc[i_min][1:]] +
            [i_b_mu[i_min](n)]
            ))
        self.reset_from_init(nrange, f_init, init)
        return nrange, np.array(init), f_init

    def needsMaxwInv(self):
        self.check()
        return any(np.diff(self._P_inv)<0)

    def reset_inv(self, iterations=100, timeout=None, kind=1, fmax=0.95,
                  npoints=None, nmax=8, adaptive=0):
        if npoints is None:
            npoints = self.npoints

        print("Resetting " + self.__repr__() + " for model " + self.Ctype.__name__.strip('_'))
        rho = []

        mu_e = []
        mu_c = []
        
        nc = []
        nlist = []
        self.warnings_inv = []
        self.status_inv = []
        if not adaptive:
            self.frange_inv = np.linspace(1e-4, fmax, npoints)
            init = arr([0.001] + [0.] * (self.n_baryon - 1) + [0.]) #All particles + mu_ch
            for i, f in enumerate(self.frange_inv):
                if sum(init) < 1e-7:
                    init = arr([1e-7] + [0.]* (self.n_baryon-1) + [0.]) #All particles + mu_ch
                # good = False
                # while not good:
                out, info = eos.stepE_rho_f(f, init, len(init) + 1, iterations, self.C, 10 + len(init)) # 10 = the info size
                fun = info[-len(init):]
                info = info[:-len(init)]
                info_dict = {'init':init, 'i': i, 'f': f, 'info': info, 'out': out, 'fun':fun}
                self.status_inv.append(info_dict)

                if int(info[6]) != 6:
                    self.warnings_inv.append(info_dict)

                if i % (npoints / 20) == 0: print('.', end='')
                rho.append(out.copy()[:-2]) ## 2 reserved for mu_c and n_c

                # rho contains all scalar fields as well
                rho[i] = np.insert(rho[i], 0, f)  #and again sp omitted
                mu_e.append(out[-2])
                mu_c.append(out[-2])
                nc.append(out[-1])
                # mu_n.append(self.mu(inv=1))
                init = out[:-1]
                _n = sum(init[:-1])
                # print(f, init, _n)
                nlist.append(_n)
                if (_n > nmax * self.n0):
                    break
        else:
            self.frange_inv = []
            f = 0.
            df = 0.01
            init = arr([0.001] + [0.00] * (self.n_baryon - 1) + [0.]) #All particles + mu_ch
            while f < fmax:
                if sum(init) < 1e-7:
                    init = arr([1e-7] + [0.]* (self.n_baryon-1) + [0.]) #All particles + mu_ch

                out = eos.stepE_rho_f(f, init, len(init) + 1, iterations, self.C)
                #len(out) = len(init) + 1 [for n_rho]

                if i % (npoints / 20) == 0: print('.', end=' ')
                rho.append(out.copy()[:-2]) ## 2 reserved for mu_c and n_c

                # rho contains all scalar fields as well
                rho[i] = np.insert(rho[i], 0, f)  #and again sp omitted
                mu_e.append(out[-2])
                mu_c.append(out[-2])
                nc.append(out[-1])
                # mu_n.append(self.mu(inv=1))
                init = out[:-1]
                _n = sum(init[:-1])
                # print(f, init, _n)
                nlist.append(_n)
                if (_n > nmax * self.n0):
                    break

# Nafiga?
        #self.frange_inv = self.frange_inv[:len(self.nrange)]
        self.frange_inv = self.frange_inv[:len(nlist)]


        #Other things left the same as for the original reset()
        self.nrange_inv = np.array(nlist)
        self.mu_c_inv = np.array(mu_c)
        self.nc_inv = np.nan_to_num(np.ascontiguousarray(arr(nc)))
        self.rho_inv = np.ascontiguousarray(arr(rho))
        self.mu_e_inv = np.ascontiguousarray(arr(mu_e))
        self._E_inv = self.Efull(rlist=self.rho_inv,
                                 nclist=self.nc_inv,
                                 mu_list=self.mu_e_inv)

        self._P_inv = np.ascontiguousarray(self.P_chem(rho=self.rho_inv,
                                                       mu_e_list=self.mu_e_inv,
                                                       nclist = self.nc_inv))

        self.set = True


    def dumpEosInv(self):
        filename = self.filenames['eos'] + '_inv'
        conc_inv = [self.rho_inv[:, i]/self.nrange_inv for i in range(1,
            self.n_baryon + 1)]
        data = np.array([
                            self.nrange_inv/self.n0,
                            self._E_inv,
                            self._P_inv,
                            self.rho_inv[:, 0]
                        ] + conc_inv + [
                            self.nc_inv / self.nrange_inv,
                            self.m_pi*(self._E_inv/self.nrange_inv - self.C.M[0]),
                            self.m_pi * self.mu_e_inv
                        ]).transpose()
        np.savetxt(join(self.foldername, filename), data, fmt='%.6f')

    def loadEosInv(self):
        filename = self.filenames['eos'] + '_inv'
        data = np.loadtxt(join(self.foldername, filename))
        self.nrange_inv = self.n0 * data[:, 0]
        self._E_inv = data[:, 1]
        self._P_inv = data[:, 2]
        self.rho_inv = data[:, 3:3 + self.n_baryon + 1]
        self.rho_inv[:, 1:] = arr([
                        r * self.nrange_inv
                        for r in self.rho_inv[:, 1:].transpose()
                    ]).transpose()
        self.nc_inv = data[:, 3+self.n_baryon + 1]
        self.mu_e_inv = data[:, 3+self.n_baryon + 2] / self.m_pi
        n_e = []
        n_mu = []
        for i, r in enumerate(self.rho_inv):
            # print(r)
            # print(eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C))
            # _mue = eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C)
            # mu_e.append(_mue)
            _mue = self.mu_e_inv[i]
            if _mue > self.m_e:
                n_e.append((_mue**2 - self.m_e**2)**(3./2) / (3 * pi **2))
            else:
                n_e.append(0.)
            if _mue > self.m_mu:
                n_mu.append((_mue**2 - self.m_mu**2)**(3./2) / (3*pi**2))
            else:
                n_mu.append(0.)

        self.n_e_inv = np.array(n_e)
        self.n_mu_inv = np.array(n_mu)

        self.set = 1

    def processMaxwInv(self, mu_init=None, show=0, shift=0, save=0):
        mu = np.nan_to_num(self.mu(inv=1)[:, 0])

        P = self._P_inv
        E = self._E_inv
        N = self.nrange_inv

        i_break = np.where(np.diff(P) < 0)[0][-1]
        i_break_mu = shift + np.where(np.diff(mu[shift:]) < 0)[0][0]
        print(i_break, i_break_mu)
        b1 = np.array([mu[:i_break_mu], P[:i_break_mu], E[:i_break_mu], N[:i_break_mu]]).transpose()
        b2 = np.array([mu[i_break+1:], P[i_break+1:], E[i_break+1:], N[i_break+1:]]).transpose()
        # kind = 'cubic'
        kind = 'linear'
        # print(b1)
        ip1 = interp1d(b1[:, 0], b1[:, 1], bounds_error=0, fill_value=0., kind=kind)
        ip2 = interp1d(b2[:, 0], b2[:, 1], bounds_error=0, fill_value=0., kind=kind)
        ie1 = interp1d(b1[:, 0], b1[:, 2], bounds_error=0, fill_value=0., kind=kind)
        ie2 = interp1d(b2[:, 0], b2[:, 2], bounds_error=0, fill_value=0., kind=kind)
        in1 = interp1d(b1[:, 0], b1[:, 3], bounds_error=0, fill_value=0., kind=kind)
        in2 = interp1d(b2[:, 0], b2[:, 3], bounds_error=0, fill_value=0., kind=kind)

        # print(b1[:i_break, 0], b2[:, 0])

        # mu_inter = np.intersect1d(b1[:, 0], b2[:, 0])
        mu_inter = np.linspace(min(b2[:, 0]), max(b1[:, 0]))
        # print(mu_inter)
        i_eq = np.argmin(abs(ip1(mu_inter) - ip2(mu_inter)))

        # print('i_eq = ', i_eq)
        if mu_init is None:
            # mu_init = [mu[i_break]]
            mu_init = [mu_inter[i_eq]]
            print('mu_init = ', mu_init)
        if show:
            plt.plot(mu_inter, ip1(mu_inter))
            plt.plot(mu_inter, ip2(mu_inter))
            plt.plot(b2[:, 0], ip2(b2[:, 0]))
            plt.plot(b1[:, 0], ip1(b1[:, 0]))
            plt.show()
        # mu_eq = bisect(lambda z: ip1(z) - ip2(z), mu_init-1., mu_init+1.)
            plt.plot(b1[:, 0], b1[:, 1])
            plt.plot(b2[:, 0], b2[:, 1])
            plt.show()
        mu_eq = root(lambda z: ip1(z) - ip2(z), mu_init, tol=1e-6).x
        print('mu_eq = ', mu_eq)
        P_eq = ip1(mu_eq)
        P_eq2 = ip2(mu_eq)
        print('P_eq = ', P_eq, 'P_eq2 = ', P_eq2)
        print('E1 = ', ie1(mu_eq))
        print('E2 = ', ie2(mu_eq))
        n1 = in1(mu_eq)
        n2 = in2(mu_eq)
        print('n1 = %.6f, n2 = %.6f' %(n1/self.n0, n2/self.n0))
        print('E1 = ', (ie1(mu_eq)/n1 - self.C.M[0]) * self.m_pi)
        print('E2 = ', (ie2(mu_eq)/n2 - self.C.M[0]) * self.m_pi)

        p1 = np.array([p for p in b1[:, 1] if p < P_eq])
        p2 = np.array([p for p in b2[:, 1] if p > P_eq])
        n1 = b1[:len(p1), 3]
        n2 = b2[(len(b2[:, 1]) - len(p2)):, 3]

        e1 = b1[:len(p1), 2]
        e2 = b2[(len(b2[:, 1]) - len(p2)):, 2]

        p_total = np.concatenate((p1, p2))
        n_total = np.concatenate((n1, n2))
        e_total = np.concatenate((e1, e2))

        self.nrange_maxw = n_total
        self._P_maxw = p_total
        self._E_maxw = e_total

        if save:
            np.savetxt(join(self.foldername, self.filenames['eos']+'_maxw'),
                    np.array([self.nrange_maxw/self.n0, self._E_maxw, self._P_maxw,
                            self.m_pi * (self._E_maxw/self.nrange_maxw
                            - self.C.M[0]) ]).transpose(),
                            fmt='%.8f')

    def check_eq(self, n, C, f, init, mu=None):
        params = eos.fun_n_eq_params()
        params.C = C
        params.n = n
        f_init = eos.dArray(1)
        f_init[0] = f
        params.f_init = f_init
        mu_min = 0.
        mu_max = 300./self.m_pi
        res = []
        for mu_e in np.linspace(mu_min, mu_max, 100):
            init[-2] = mu_e
            res.append(eos.func_f_eq_rho_anal())

    def Efull(self, n=None, nc=None, mu_e=None, rlist=None, mu_list=None, nclist=None):
        if (n is None) or (mu_e is None) or (nc is None):
            if rlist is None:
                rlist = self.rho
            if mu_list is None:
                mu_list = self.mu_e
            if nclist is None:
                nclist = self.nc
            _E = []
            for i, n in enumerate(rlist):
                mu_e = mu_list[i]
                E = eos.E_rho(n, mu_e, self.C)
                n_e = 0.
                n_mu = 0.
                if (mu_e ** 2 - self.m_e ** 2 > 0):
                    n_e += (mu_e ** 2 - self.m_e ** 2) ** (1.5) / (3 * pi ** 2)

                if (mu_e ** 2 - self.m_mu ** 2 > 0):
                    n_mu += (mu_e ** 2 - self.m_mu ** 2) ** (1.5) / (3 * pi ** 2)

                E += eos.kineticInt(n_e, self.m_e, 2.)
                E += eos.kineticInt(n_mu, self.m_mu, 2.)
                E += nclist[i] * mu_e
                _E.append(E)

            return np.array(_E)
        else:
            epart = np.zeros((9), dtype='float')

            E = eos.E_rho(n, mu_e, self.C)
            n_e = 0.
            n_mu = 0.
            if (mu_e ** 2 - self.m_e ** 2 > 0):
                n_e += (mu_e ** 2 - self.m_e ** 2) ** (1.5) / (3 * pi ** 2)

            if (mu_e ** 2 - self.m_mu ** 2 > 0):
                n_mu += (mu_e ** 2 - self.m_mu ** 2) ** (1.5) / (3 * pi ** 2)

            E += nc * mu_e
            E += eos.kineticInt(n_e, self.m_e, 2.)
            E += eos.kineticInt(n_mu, self.m_mu, 2.)

            return E, epart


    def Efull2(self, n=None, mu_e=None):
        if (n is None) or (mu_e is None):
            _E = []
            for i, n in enumerate(self.rho):
                mu_e = self.mu_e[i]
                E = eos.E_rho(n, mu_e, self.C)
                n_e = 0.
                n_mu = 0.
                if (mu_e ** 2 - self.m_e ** 2 > 0):
                    n_e += (mu_e ** 2 - self.m_e ** 2) ** (1.5) / (3 * pi ** 2)

                if (mu_e ** 2 - self.m_mu ** 2 > 0):
                    n_mu += (mu_e ** 2 - self.m_mu ** 2) ** (1.5) / (3 * pi ** 2)

                E += eos.kineticInt(n_e, self.m_e, 2.)
                E += eos.kineticInt(n_mu, self.m_mu, 2.)
                _E.append(E)

            return np.array(_E)
        else:
            epart = np.zeros((9), dtype='float')

            E = eos.E_rho(n, mu_e, self.C)
            n_e = 0.
            n_mu = 0.
            if (mu_e ** 2 - self.m_e ** 2 > 0):
                n_e += (mu_e ** 2 - self.m_e ** 2) ** (1.5) / (3 * pi ** 2)

            if (mu_e ** 2 - self.m_mu ** 2 > 0):
                n_mu += (mu_e ** 2 - self.m_mu ** 2) ** (1.5) / (3 * pi ** 2)

            E += eos.kineticInt(n_e, self.m_e, 2.)
            E += eos.kineticInt(n_mu, self.m_mu, 2.)

            return E, epart

    def P_chem2(self, leptons=True):
        res = []
        sp = 1
        nlist = self.rho
        for i, n in enumerate(nlist):
            n = arr(n)
            mu_e = self.mu_e[i]
            n_e = 0.
            n_mu = 0.
            if (mu_e ** 2 - self.m_e ** 2 > 0):
                n_e += (mu_e ** 2 - self.m_e ** 2) ** (1.5) / (3 * pi ** 2)

            if (mu_e ** 2 - self.m_mu ** 2 > 0):
                n_mu += (mu_e ** 2 - self.m_mu ** 2) ** (1.5) / (3 * pi ** 2)

            if leptons:
                E, epart = self.Efull2(n=n, mu_e=self.mu_e[i])

            else:
                E = eos.E_rho(n, self.mu_e[i], self.C)

            sum = 0.

            for j, r in enumerate(n[sp:]):
                sum += r * eos.mu_rho(n, j + sp, mu_e, self.C)

            if leptons:
                sum += n_e * mu_e + n_mu * mu_e

            res.append(sum - E)
        return np.array(res) * self.mpi4_2_mevfm3


    def P_chem(self, rho=None, mu_e_list=None, leptons=True, nclist=None):
        res = []
        sp = 1
        if rho is None:
            nlist = self.rho
        else:
            nlist = rho

        if nclist is None:
            nclist = self.nc

        if mu_e_list is None:
            mu_e_list = self.mu_e

        for i, n in enumerate(nlist):
            n = arr(n)
            mu_e = mu_e_list[i]
            n_e = 0.
            n_mu = 0.
            if (mu_e ** 2 - self.m_e ** 2 > 0):
                n_e += (mu_e ** 2 - self.m_e ** 2) ** (1.5) / (3 * pi ** 2)

            if (mu_e ** 2 - self.m_mu ** 2 > 0):
                n_mu += (mu_e ** 2 - self.m_mu ** 2) ** (1.5) / (3 * pi ** 2)

            if leptons:
                E, epart = self.Efull(n=n, nc=nclist[i], mu_e=mu_e_list[i])
                # E += self.nc[i] * mu_e
            else:
                E = eos.E_rho(n, self.mu_e[i], self.C)

            # sum = 0.
            sum = mu_e * nclist[i]

            for j, r in enumerate(n[sp:]):
                sum += r * eos.mu_rho(n, j + sp, mu_e, self.C)

            if leptons:
                sum += n_e * mu_e + n_mu * mu_e

            res.append(sum - E)
        return np.array(res) * self.mpi4_2_mevfm3


    def get_upper_branch(self):
        self.loadEos()
        out = self.checkMultiSol()
        if out == []:
            print("No other branches found")
            return
        f_init = np.array([out[-1][-1]])
        i_start = len(self.nrange) - 1
        # Third solution with selfconsistent concentrations usually appears
        # earlier than it was with "old" concentrations => factor 0.75
        i_end = int(0.75 * out[1][1])
        print(i_end)
        n_tosolve = self.nrange[i_start:i_end:-1]
        n_init = self.rho[i_start]
        init = np.insert(n_init[2:], -1, self.mu_e[i_start])

        # Solving and populating concentrations
        flist = []
        rlist = []
        mu_e = []
        self.nc2 = []
        for i_solve, _n in enumerate(n_tosolve):
        #     print(_n/m.n0)
            res = eos.stepE_rho(_n, init, f_init, len(init), 30, self.mu_e[i_solve], self.C)
            n_S = np.array([_n - res[0], res[0]])
            self.nc2.append(res[-1])
            f = eos.f_eq_rho(n_S, f_init, 1, res[-2], self.C)
        #     f = eos.f_eq(n_S, f_init, 1, m.C)
        #     print(res, f)
            rlist.append([_n - res[0], res[0]])
            flist.append(f[0])
            f_init = f
            mu_e.append(res[-2])
        mu_e =np.array(mu_e)
        rlist = np.array(rlist)
        flist = np.array(flist)
        self.nc2 = np.array(self.nc2)
        self.rho_upper = np.insert(rlist, 0, flist.transpose(), axis=1)
        self.mu_e_upper = mu_e
        self.n_upper = n_tosolve

        # print(self..checkMultiSol)

        # Stripping the new concentrations for the actual solution

    def refine_upper(self):
        sols = self.checkMultiSol(nrange=self.n_upper, rho=self.rho_upper,
            mu=self.mu_e_upper, npoints=len(self.n_upper), ret_3br=0)
        # print(sols)
        i_break = sols[-1][1]
        print(i_break)
        self.rho_upper2 = self.rho_upper[:i_break]
        self.mu_e_upper2 = self.mu_e_upper[:i_break]
        self.n_upper2 = self.n_upper[:i_break]
        self.nc2 = self.nc2[:i_break]
        # self.E_upper = np.array([eos.E_rho(n, self.mu_e_upper2[i], self.C) for i, n in enumerate(self.rho_upper2)])
        self._E2 = self._E.copy()
        self.rho2 = self.rho.copy()
        self.mu_e2 = self.mu_e.copy()
        self.E_upper = self.Efull(rlist=self.rho_upper2,
            mu_list=self.mu_e_upper2, nclist=self.nc2)
        self._P2 = self._P.copy()
        self.nclist2 = self.nc.copy()
        self.P_upper = self.P_chem(rho=self.rho_upper2, mu_e_list=self.mu_e_upper2, nclist=self.nc2)
        shift = len(self._E2) - len(self.E_upper)
        flag = 0 #VERY DIRTY!!!!!!!!!!!!!!!!!!!! assumes that third solution
                                                #doesn't disappear and is energetically
                                                #favorable for all densities!
        for i, e in enumerate(self.E_upper[::-1]):
            if self._E2[shift + i] > e:
                flag = 1
            if flag:
                self._E2[shift + i] = e
                self.rho2[shift + i] = self.rho_upper2[::-1][i]
                self.mu_e2[shift + i] = self.mu_e_upper2[::-1][i]
                self._P2[shift + i] = self.P_upper[::-1][i]
                self.nclist2[shift + i] = self.nc2[::-1][i]

class Rcc(Rho):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.foldername = self.foldername + '_chi'
        self.basefoldername = self.foldername
        if not os.path.exists(self.foldername):
            os.makedirs(self.foldername)
        self.C.chi_r_prime = 1

class Rcp(Rho):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.foldername = self.foldername + '_phi'
        self.basefoldername = self.foldername
        if not os.path.exists(self.foldername):
            os.makedirs(self.foldername)
        self.C.chi_r_prime = 2


class RhoNucleon(Rho, Nucleon):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['rcond'] = 'rcond_nucl.dat'
        self.filenames['mass_crust'] = 'rc_masses_N.dat'
        self.filenames['n_rho'] = 'nr_nucl.dat'


class RhoHyper(Rho, Hyperon):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['rcond'] = 'rcond_hyper.dat'
        self.filenames['mass_crust'] = 'rc_masses_H.dat'
        self.filenames['n_rho'] = 'nr_hyper.dat'


class RhoHyperPhi(Rho, HyperonPhi):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['rcond'] = 'rcond_hyper_phi.dat'
        self.filenames['mass_crust'] = 'rc_masses_Hp.dat'
        self.filenames['n_rho'] = 'nr_hyper_phi.dat'


class RhoHyperPhiSigma(Rho, HyperonPhiSigma):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['rcond'] = 'rcond_hyper_phi_sigma.dat'
        self.filenames['mass_crust'] = 'rc_masses_Hps.dat'
        self.filenames['n_rho'] = 'nr_hyper_phi_sigma.dat'


class RccNucleon(Rcc, Nucleon):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['rcond'] = 'rcond_nucl.dat'
        self.filenames['mass_crust'] = 'rc_masses_N.dat'
        self.filenames['n_rho'] = 'nr_nucl.dat'


class RccHyper(Rcc, Hyperon):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['rcond'] = 'rcond_hyper.dat'
        self.filenames['mass_crust'] = 'rc_masses_H.dat'
        self.filenames['n_rho'] = 'nr_hyper.dat'

class RccHyperPhi(Rcc, HyperonPhi):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['rcond'] = 'rcond_hyper_phi.dat'
        self.filenames['mass_crust'] = 'rc_masses_Hp.dat'
        self.filenames['n_rho'] = 'nr_hyper_phi.dat'


class RccHyperPhiSigma(Rcc, HyperonPhiSigma):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['rcond'] = 'rcond_hyper_phi_sigma.dat'
        self.filenames['mass_crust'] = 'rc_masses_Hps.dat'
        self.filenames['n_rho'] = 'nr_hyper_phi_sigma.dat'

class RcpNucleon(Rcp, Nucleon):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['rcond'] = 'rcond_nucl.dat'
        self.filenames['mass_crust'] = 'rc_masses_N.dat'


class RcpHyper(Rcp, Hyperon):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['rcond'] = 'rcond_hyper.dat'
        self.filenames['mass_crust'] = 'rc_masses_H.dat'


class RcpHyperPhi(Rcp, HyperonPhi):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['rcond'] = 'rcond_hyper_phi.dat'
        self.filenames['mass_crust'] = 'rc_masses_Hp.dat'


class RcpHyperPhiSigma(Rcp, HyperonPhiSigma):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['rcond'] = 'rcond_hyper_phi_sigma.dat'
        self.filenames['mass_crust'] = 'rc_masses_Hps.dat'





class DeltaOnly(DeltaBase, Hyperon):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.foldername = join(self.foldername, 'DOnly')
        if hasattr(self.C, 'setDeltaConstants'):
            self.C.setDeltaOnlyConstants()
            self.n_baryon = 6
            self.part_names = ['n', 'p', 'D-', 'D0', 'D+', 'D++']
            self.filenames['eos'] = 'eos_donly.dat'
            self.filenames['meff'] = 'meff_donly.dat'
            self.filenames['vs'] = 'vs_donly.dat'
            self.filenames['grig'] = 'grig_donly.dat'
            self.filenames['mu'] = 'mu_donly.dat'
            self.filenames['fortin'] = 'fortin_donly.dat'
            self.filenames['pf'] = 'pf_donly.dat'
            self.filenames['mass_crust'] = 'mass_crust_donly.dat'
            self.filenames['mass_nocrust'] = 'mass_nocrust_donly.dat'
            self.filenames['S'] = 'S_donly.dat'
            self.filenames['V'] = 'V_donly.dat'
            self.filenames['I'] = 'I_donly.dat'
            self.filenames['density'] = 'dens_donly.dat'
            self.filenames['eta'] = 'eta_donly.dat'
            self.filenames['var_u'] = 'nc_do_U.dat'
            self.filenames['var_xo'] = 'nc_do_xo.dat'
            self.filenames['var_xr'] = 'nc_do_xr.dat'


    def loadEos(self):
        try:
            eos_tab = np.loadtxt(join(self.foldername, self.filenames['eos']),
                skiprows=1)
            self.nrange = eos_tab[:, 0] * self.n0
            self._E = eos_tab[:, 1]
            self._P = eos_tab[:, 2]
            self.set = 1
            rhos = self.nrange*eos_tab[:, 3: 3 + self.n_baryon].transpose()
            # print(rhos)
            self.rho = np.insert(rhos.transpose(), 0, eos_tab[:, 3 + self.n_baryon], axis=1)
            mu_e = []
            n_e = []
            n_mu = []
            for r in self.rho:
                # print(r)
                # print(eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C))
                _mue = eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C)
                mu_e.append(_mue)
                if _mue > self.m_e:
                    n_e.append((_mue**2 - self.m_e**2)**(3./2) / (3 * pi **2))
                else:
                    n_e.append(0.)
                if _mue > self.m_mu:
                    n_mu.append((_mue**2 - self.m_mu**2)**(3./2) / (3*pi**2))
                else:
                    n_mu.append(0.)

            self.n_e = np.array(n_e)
            self.n_mu = np.array(n_mu)

            # print(self.n_e +self.n_mu - self.rho[:,2])
            # exit()

            self.mu_e = np.array(mu_e)

            self.xDUp = ((3*pi**2 * self.rho[:, 1])**(1./3) - (3*pi**2 * self.n_e)**(1./3))**(3.) / (3*pi**2 * self.nrange)
        except FileNotFoundError:
            raise
            # self.dumpEos()


class DeltaSym2(DeltaOnly):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['eos'] = 'eos_dsym2.dat'
        self.filenames['meff'] = 'meff_dsym2.dat'
        self.filenames['vs'] = 'vs_dsym2.dat'
        self.filenames['grig'] = 'grig_dsym2.dat'
        self.filenames['mu'] = 'mu_dsym2.dat'
        self.filenames['fortin'] = 'fortin_dsym2.dat'
        self.filenames['pf'] = 'pf_dsym2.dat'
        self.filenames['mass_crust'] = 'mass_crust_dsym2.dat'
        self.filenames['mass_nocrust'] = 'mass_nocrust_dsym2.dat'
        self.filenames['S'] = 'S_dsym2.dat'
        self.filenames['V'] = 'V_dsym2.dat'
        self.filenames['I'] = 'I_dsym2.dat'
        self.filenames['density'] = 'dens_dsym2.dat'
        self.filenames['eta'] = 'eta_dsym2.dat'
        # self.filenames['var_u'] = 'nc_do_U.dat'
        # self.filenames['var_xo'] = 'nc_do_xo.dat'
        # self.filenames['var_xr'] = 'nc_do_xr.dat'


    def checkDeltaMultiSol(self, nrange=None, npoints=20, break_lost=0):
        # if mu is None:
        #     mu = self.mu
        if nrange is None:
            nrange = self.nrange
        nd_res = []
        out = []
        i_check = range(0, len(nrange), len(nrange)//npoints)
        for i in i_check:
            nd, res = self.get_fn(nrange[i], f_init=0.5)
            # n_out = [] #densities at which more that 1 solution occurs
            # f_out = [] #approximate values of the roots for these densities
            sch_prev = 0
            sign_changes = 0
            r_prev = res[1]
            for k, r in enumerate(res):
                if r*r_prev < 0:
                    sign_changes += 1
                    print(i, sign_changes, nd[k])
                    nd_res.append(nd[k])
                    r_prev = r

            if sign_changes > 1:
                out.append([sign_changes, i, nrange[i]] + nd_res)
                return out

            if sign_changes < sch_prev:
                if break_lost:
                    return out
            sch_prev = sign_changes
        return out

    def reset(self, iterations=30, timeout=None):

        print("Resetting " + self.__repr__() + " for model " + self.Ctype.__name__.strip('_'))
        init = arr([0.])

        rho = []

        mu_e = []
        f = arr([0.])  # sp omitted
        for i, _n in enumerate(self.nrange):

            if timeout is None:
                init = eos.stepE_dsym(_n, init, f, len(init), iterations, self.C)
            else:
                queue = Queue()
                p = Process(target=self.stepE_dsym, args=(_n, init, f, len(init),
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
                    self.set = 1
                    return
                init = queue.get(timeout=None)
                # print init

            if i % (len(self.nrange) / 20) == 0:
                print('.', end=' ')
            n_d = init[0]
            n_n = (_n - n_d)/2
            _rho = arr([n_n, n_n, n_d/4, n_d/4, n_d/4, n_d/4])
            rho.append(np.ascontiguousarray(_rho))
            # rho[i] = np.insert(rho[i], 0, _n - np.sum(init))
            f = eos.f_eq(rho[i], f, 1, self.C)  # sp omitted
            if self.verbose:
                pass  # TODO: some verbose output
            # rho contains all scalar fields as well
            rho[i] = np.insert(rho[i], 0, f)  # and again sp omitted
            mu_e.append(eos.mu(rho[i], 1, self.C) -
                        eos.mu(rho[i], 2, self.C))

        self.rho = np.ascontiguousarray(arr(rho))
        self.mu_e = np.ascontiguousarray(arr(mu_e))
        eparts = []
        _E = []
        for z in self.rho:
            epart = np.zeros((9), dtype='float')
            _E.append(eos.E(z, self.C))
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

    def get_fn(self, n, f_init=0.5):
        nlist = np.linspace(0, n, 500)
        res = []
        finit = f_init
        for _n in nlist:
            eq = eos.wrap_fun_dsym(n, _n, finit, self.C)
            res.append(eq)
        return np.array(nlist), np.array(res)

    def get_fn_nd(self, nd, nlist, f_init=0.5):
        res = []
        finit = f_init
        for _n in nlist:
            eq = eos.wrap_fun_dsym(_n, nd, finit, self.C)
            res.append(eq)
        return np.array(res)

    def get_fn_f(self, n, f):
        nlist = np.linspace(0, n, 100)
        res = []
        finit = f
        for _n in nlist:
            eq = eos.wrap_fun_dsym_f(n, _n, finit, self.C)
            res.append(eq)
        return np.array(nlist), np.array(res)

    def inspect_eq(self):
        self.check()


        fig, ax = plt.subplots()

        # func = list(map(
        #     lambda z: func(z, swArr(self.rho[0, 1:])), frange
        # ))

        n = max(self.nrange)
        res0 = self.get_fn(0.01)
        l, = ax.plot(res0[0]/0.01, res0[1])

        plt.subplots_adjust(left=0.25, bottom=0.25)
        axn = plt.axes([0.25, 0.1, 0.65, 0.03])
        sn = Slider(axn, 'N', 0, len(self.nrange), valinit=1.)

        def update(val):
            print('val=', val)
            _n = self.nrange[int(val)]
            res = self.get_fn(_n, f_init=self.rho[int(val), 0])
            l.set_ydata(res[1])
            l.set_xdata(res[0]/_n)
            fig.canvas.draw_idle()

        sn.on_changed(update)
        ax.set_ylim([-5, 5])
        plt.show()

    def get_proper_second_fucking_branch(self, i_start=0, init=0.5):
        init = np.array([init])
        f = np.array([self.rho[i_start, 0]])
        rho = []
        for n in self.nrange[i_start:]:
        #     _res = m.getSymm(n, lastx=lastx, lastf=lastf)
            _init = eos.stepE_dsym(n, init, f, 1, 30, self.C)
            fun = eos.wrap_fun_dsym(n, _init[0], f[0], self.C)
            if (abs(fun) < 1e-4) and (abs(fun) > 0.0):
                init = _init
            else:
                _init[0] = 0.

            n_d = _init[0]
            n_n = (n - n_d)/2
            _rho = np.array([n_n, n_n, n_d/4, n_d/4, n_d/4, n_d/4])
            # rho[i] = np.insert(rho[i], 0, _n - np.sum(init))
            f = eos.f_eq(_rho, f, 1, self.C)
            _rho = np.insert(_rho, 0, f)
            rho.append(np.ascontiguousarray(_rho))
            print(n, init[0], n_d, f[0], fun)
        self.rho = np.nan_to_num(np.ascontiguousarray(rho))
        self._E = np.nan_to_num(np.array([eos.E(z, self.C) for z in self.rho]))
        self._P = np.nan_to_num(self.P_chem(rho))
        self.set = 1


    def setSecondBranch():
        pass



class DeltaPhi(DeltaBase, HyperonPhi):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['var_u'] = 'nc_dp_U.dat'
        self.filenames['var_xo'] = 'nc_dp_xo.dat'
        self.filenames['var_xr'] = 'nc_dp_xr.dat'

class DeltaPhi2(DeltaBase, HyperonPhi2):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)


class DeltaPhiSigma(DeltaBase, HyperonPhiSigma):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['var_u'] = 'nc_dps_U.dat'
        self.filenames['var_xo'] = 'nc_dps_xo.dat'
        self.filenames['var_xr'] = 'nc_dps_xr.dat'

class RhoDeltaPhi(DeltaPhi, Rho):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.basefoldername = join(self.basefoldername, 'rcond')
        self.filenames['rcond'] = 'rcond_hyper_phi.dat'
        self.filenames['mass_crust'] = 'rc_masses_Hp.dat'
        self.filenames['n_rho'] = 'nr_dp.dat'


    def mu(self, nrange=None, branch_3=0, inv=0):
        if nrange is None:
            nrange = self.nrange

        if not inv:
            self.check()

        if inv:
            conc = self.rho_inv
            mu_e = self.mu_e_inv
        else:
            if branch_3:
                conc = self.rho2
                mu_e = self.mu_e
            else:
                conc = self.rho
                mu_e = self.mu_e


        mu = []
        for j, n in enumerate(conc):
            mu.append([eos.mu_rho(n, i+1, mu_e[j], self.C) for i in range(self.n_baryon)])
        return arr(mu)


    def dumpEos(self):
        super().dumpEos()

        rcond = arr([self.nrange/self.n0, self.nc/self.nrange]).transpose()
        np.savetxt(join(self.foldername, self.filenames['rcond'])
                   ,rcond, fmt='%.8f')
        

    def checkEq(self):
        params = eos.fun_n_eq_params()
        params.C = self.C
        params.dimF_init = 1
        params.misc = 0.
        for i, r in enumerate(self.rho):
            params

    def loadEos(self):
        super().loadEos()

        rcond = np.loadtxt(join(self.foldername, self.filenames['rcond']))
        self.nc = self.nrange * rcond[:, 1]

        mu_e = []
        n_e = []
        n_mu = []

        grig = np.loadtxt(join(self.foldername, self.filenames['grig']),
            skiprows=1)
        mu_e = grig[:, 13]/self.m_pi
        # print(mu_e)
        # exit()

        n_e = []
        n_mu = []
        for i, r in enumerate(self.rho):
            # print(r)
            # print(eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C))
            # _mue = eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C)
            # mu_e.append(_mue)
            _mue = mu_e[i]
            if _mue > self.m_e:
                n_e.append((_mue**2 - self.m_e**2)**(3./2) / (3 * pi **2))
            else:
                n_e.append(0.)
            if _mue > self.m_mu:
                n_mu.append((_mue**2 - self.m_mu**2)**(3./2) / (3*pi**2))
            else:
                n_mu.append(0.)

        self.n_e = np.array(n_e)
        self.n_mu = np.array(n_mu)

        # print(self.n_e +self.n_mu - self.rho[:,2])
        # exit()

        self.mu_e = np.array(mu_e)
        self.mu_c = np.array(mu_e)
        self.n_e = np.array(n_e)
        self.n_mu = np.array(n_mu)

        # print(self.n_e +self.n_mu - self.rho[:,2])
        # exit()

        self.mu_e = np.array(mu_e)
        self.mu_c = np.array(mu_e)


class RhoDeltaOnly(DeltaOnly, Rho):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['rcond'] = 'rcond_donly.dat'
        self.filenames['mass_crust'] = 'rc_masses_donly.dat'
        self.filenames['n_rho'] = 'nr_do.dat'

    def mu(self, nrange=None, branch_3=0, inv=0):
        if nrange is None:
            nrange = self.nrange

        if not inv:
            self.check()

        if inv:
            conc = self.rho_inv
            mu_e = self.mu_e_inv
        else:
            if branch_3:
                conc = self.rho2
                mu_e = self.mu_e
            else:
                conc = self.rho
                mu_e = self.mu_e

        mu = []
        for j, n in enumerate(conc):
            mu.append([eos.mu_rho(n, i+1, mu_e[j], self.C) for i in range(self.n_baryon)])
        return arr(mu)


    def dumpEos(self):
        super().dumpEos()

        rcond = arr([self.nrange/self.n0, self.nc/self.nrange]).transpose()
        np.savetxt(join(self.foldername, self.filenames['rcond'])
                   ,rcond, fmt='%.8f')

    def loadEos(self):
        super().loadEos()

        rcond = np.loadtxt(join(self.foldername, self.filenames['rcond']))
        self.nc = self.nrange * rcond[:, 1]

        mu_e = []
        n_e = []
        n_mu = []

        grig = np.loadtxt(join(self.foldername, self.filenames['grig']),
            skiprows=1)
        mu_e = grig[:, 13]/self.m_pi
        # print(mu_e)
        # exit()

        n_e = []
        n_mu = []
        for i, r in enumerate(self.rho):
            # print(r)
            # print(eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C))
            # _mue = eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C)
            # mu_e.append(_mue)
            _mue = mu_e[i]
            if _mue > self.m_e:
                n_e.append((_mue**2 - self.m_e**2)**(3./2) / (3 * pi **2))
            else:
                n_e.append(0.)
            if _mue > self.m_mu:
                n_mu.append((_mue**2 - self.m_mu**2)**(3./2) / (3*pi**2))
            else:
                n_mu.append(0.)

        self.n_e = np.array(n_e)
        self.n_mu = np.array(n_mu)

        # print(self.n_e +self.n_mu - self.rho[:,2])
        # exit()

        self.mu_e = np.array(mu_e)
        self.mu_c = np.array(mu_e)
        self.n_e = np.array(n_e)
        self.n_mu = np.array(n_mu)

        # print(self.n_e +self.n_mu - self.rho[:,2])
        # exit()

        self.mu_e = np.array(mu_e)
        self.mu_c = np.array(mu_e)


class RhoDeltaPhiSigma(DeltaPhiSigma, Rho):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.basefoldername = join(self.basefoldername, 'rcond')
        self.filenames['rcond'] = 'rcond_hyper_phi_sigma.dat'
        self.filenames['mass_crust'] = 'rc_masses_Hps.dat'
        self.filenames['n_rho'] = 'nr_dps.dat'

    def mu(self, nrange=None, branch_3=0, inv=0):
        if nrange is None:
            nrange = self.nrange

        if not inv:
            self.check()

        if inv:
            conc = self.rho_inv
            mu_e = self.mu_e_inv
        else:
            if branch_3:
                conc = self.rho2
                mu_e = self.mu_e
            else:
                conc = self.rho
                mu_e = self.mu_e


        mu = []
        for j, n in enumerate(conc):
            mu.append([eos.mu_rho(n, i+1, mu_e[j], self.C) for i in range(self.n_baryon)])
        return arr(mu)


    def dumpEos(self):
        super().dumpEos()

        rcond = arr([self.nrange/self.n0, self.nc/self.nrange]).transpose()
        np.savetxt(join(self.foldername, self.filenames['rcond'])
                   ,rcond, fmt='%.8f')

    def loadEos(self):
        super().loadEos()

        rcond = np.loadtxt(join(self.foldername, self.filenames['rcond']))
        self.nc = self.nrange * rcond[:, 1]

        mu_e = []
        n_e = []
        n_mu = []

        grig = np.loadtxt(join(self.foldername, self.filenames['grig']),
            skiprows=1)
        mu_e = grig[:, 13]/self.m_pi
        # print(mu_e)
        # exit()

        n_e = []
        n_mu = []
        for i, r in enumerate(self.rho):
            # print(r)
            # print(eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C))
            # _mue = eos.mu(r, 1, self.C) - eos.mu(r, 2, self.C)
            # mu_e.append(_mue)
            _mue = mu_e[i]
            if _mue > self.m_e:
                n_e.append((_mue**2 - self.m_e**2)**(3./2) / (3 * pi **2))
            else:
                n_e.append(0.)
            if _mue > self.m_mu:
                n_mu.append((_mue**2 - self.m_mu**2)**(3./2) / (3*pi**2))
            else:
                n_mu.append(0.)

        self.n_e = np.array(n_e)
        self.n_mu = np.array(n_mu)

        # print(self.n_e +self.n_mu - self.rho[:,2])
        # exit()

        self.mu_e = np.array(mu_e)
        self.mu_c = np.array(mu_e)
        self.n_e = np.array(n_e)
        self.n_mu = np.array(n_mu)

        # print(self.n_e +self.n_mu - self.rho[:,2])
        # exit()

        self.mu_e = np.array(mu_e)
        self.mu_c = np.array(mu_e)



class RccDeltaOnly(Rcc, DeltaOnly):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['rcond'] = 'rcond_donly.dat'
        self.filenames['mass_crust'] = 'rc_masses_donly.dat'
        self.filenames['n_rho'] = 'nr_do.dat'

class RccDeltaPhi(Rcc, DeltaPhi):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['rcond'] = 'rcond_hyper_phi.dat'
        self.filenames['mass_crust'] = 'rc_masses_Hp.dat'
        self.filenames['n_rho'] = 'nr_dp.dat'

class RccDeltaPhiSigma(Rcc, DeltaPhiSigma):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['rcond'] = 'rcond_hyper_phi_sigma.dat'
        self.filenames['mass_crust'] = 'rc_masses_Hps.dat'
        self.filenames['n_rho'] = 'nr_dps.dat'

class RcpDeltaPhi(Rcp, DeltaPhi):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['rcond'] = 'rcond_hyper_phi.dat'
        self.filenames['mass_crust'] = 'rc_masses_Hp.dat'

class RcpDeltaPhiSigma(Rcp, DeltaPhiSigma):
    def __init__(self, C, basefolder_suffix=''):
        super().__init__(C, basefolder_suffix=basefolder_suffix)
        self.filenames['rcond'] = 'rcond_hyper_phi_sigma.dat'
        self.filenames['mass_crust'] = 'rc_masses_Hps.dat'



class Model(Wrapper):
    def __init__(self, C=eos.KVOR, K0=None, f0=None, J0=None, suffix=None, basefolder_suffix=''):
        if any([K0, f0, J0]):
            params = []
            names = []
            #a dirty workaround
            for p in [['K0', K0], ['J0', J0], ['f0', f0]]:
                if p[1] is not None:
                    names.append(p[0])
                    params.append(p[1])
                    basefolder_suffix+=p[0]+"%.2f"%(p[1])

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
        self.delta_phi2 = DeltaPhi2(C, basefolder_suffix=basefolder_suffix)
        self.delta_sym = DeltaSym(C, basefolder_suffix=basefolder_suffix)
        self.delta_asym = DeltaAsym(C, basefolder_suffix=basefolder_suffix)
        self.delta_sym2 = DeltaSym2(C, basefolder_suffix=basefolder_suffix)
        self.delta_only = DeltaOnly(C, basefolder_suffix=basefolder_suffix)
        self.children = [self.nucl, self.sym, self.neutr, self.hyper,
                         self.hyper_phi, self.hyper_phi_sigma,
                         self.delta, self.delta_phi, self.delta_phi2, self.delta_sym, self.delta_only]

        self.rcond_nucl = RhoNucleon(C, basefolder_suffix=basefolder_suffix)
        self.rcond_hyper = RhoHyper(C, basefolder_suffix=basefolder_suffix)
        self.rcond_hyper_phi = RhoHyperPhi(C, basefolder_suffix=basefolder_suffix)
        self.rcond_hyper_phi_sigma = RhoHyperPhiSigma(C, basefolder_suffix=basefolder_suffix)
        self.rcond_delta_phi = RhoDeltaPhi(C, basefolder_suffix=basefolder_suffix)
        self.rcond_delta_phi_sigma = RhoDeltaPhiSigma(C, basefolder_suffix=basefolder_suffix)
        self.rcond_delta_only = RhoDeltaOnly(C, basefolder_suffix=basefolder_suffix)

        self.rcc_nucl = RccNucleon(C, basefolder_suffix=basefolder_suffix)
        self.rcc_hyper = RccHyper(C, basefolder_suffix=basefolder_suffix)
        self.rcc_hyper_phi = RccHyperPhi(C, basefolder_suffix=basefolder_suffix)
        self.rcc_hyper_phi_sigma = RccHyperPhiSigma(C, basefolder_suffix=basefolder_suffix)
        self.rcc_delta_only = RccDeltaOnly(C, basefolder_suffix=basefolder_suffix)
        self.rcc_delta_phi = RccDeltaPhi(C, basefolder_suffix=basefolder_suffix)
        self.rcc_delta_phi_sigma = RccDeltaPhiSigma(C, basefolder_suffix=basefolder_suffix)

        self.rcp_nucl = RcpNucleon(C, basefolder_suffix=basefolder_suffix)
        self.rcp_hyper = RcpHyper(C, basefolder_suffix=basefolder_suffix)
        self.rcp_hyper_phi = RcpHyperPhi(C, basefolder_suffix=basefolder_suffix)
        self.rcp_hyper_phi_sigma = RcpHyperPhiSigma(C, basefolder_suffix=basefolder_suffix)
        self.rcp_delta_phi = RcpDeltaPhi(C, basefolder_suffix=basefolder_suffix)
        self.rcp_delta_phi_sigma = RcpDeltaPhiSigma(C, basefolder_suffix=basefolder_suffix)

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
        #     if not os.path.exists(self.foldername):   o
        #         os.makedirs(self.foldername)
        #     for child in self.children:
        #         child.foldername = self.foldername
        self.filenames['eos'] = 'eos.dat'
        self.filenames.update(eta='eta_F.dat')
        self.filenames.update(props='props.dat')
        self.filenames.update(intro='intro.dat')

    def setDeltaPotential(self, U):
        xs_d = self.getDeltaXs(U)
        xo = self.delta_phi.C.X_o[10]
        xr = self.delta_phi.C.X_r[10]
        self.setDeltaConst(np.array([xs_d for i in range(4)]),
             np.array([xo for i in range(4)]),
             np.array([xr for i in range(4)]),
             'xs=%.2f U = %.2f'%(xs_d, U))

    def setDeltaXo(self, xo):
        xs = self.delta_phi.C.X_s[10]
        xr = self.delta_phi.C.X_r[10]
        S, V = self.dumpPotentials()
        m = self.delta
        iS = interp1d(m.nrange/m.n0, S)
        iV = interp1d(m.nrange/m.n0, V)
        U = iS(1.) * xs + iV(1.) * xo
        self.setDeltaConst(np.array([xs for i in range(4)]),
             np.array([xo for i in range(4)]),
             np.array([xr for i in range(4)]),
             'xo=%.2f U = %.2f'%(xo, U))

    def setDeltaXr(self, xr):
        xs = self.delta_phi.C.X_s[10]
        xo = self.delta_phi.C.X_o[10]
        self.setDeltaConst(np.array([xs for i in range(4)]),
             np.array([xo for i in range(4)]),
             np.array([xr for i in range(4)]),
             'xr=%.2f'%(xr))

    def getDeltaXs(self, U):
        i = 8
        S, V = self.getPots(self.n0)
        m = self.delta
        # iS = interp1d(m.nrange/m.n0, np.nan_to_num(S), kind='linear')
        # iV = interp1d(m.nrange/m.n0, np.nan_to_num(V), kind='linear')
        xs_d = (U - self.delta_phi.C.X_o[i]*V)/S

        return xs_d

    def setDeltaConst(self, Xs, Xo, Xr, folder_suffix):
        for m in [self.delta, self.delta_phi,self.delta_phi_sigma,
                  self.delta_phi2, self.delta_sym, self.delta_sym2, self.delta_only,
                  self.delta_asym,
                  self.rcond_delta_phi, self.rcond_delta_phi_sigma,
                  self.rcond_delta_only,
                  self.rcc_delta_only,
                  self.rcc_delta_phi, self.rcc_delta_phi_sigma,
                  self.rcp_delta_phi, self.rcp_delta_phi_sigma]:
            m.C.setDeltaRho(Xr)
            m.C.setDeltaOmega(Xo)
            m.C.setDeltaSigma(Xs)
            m.foldername = m.basefoldername + folder_suffix
            if not os.path.exists(m.foldername):
                os.makedirs(m.foldername)
            m.set = 0

    def setParams(self, Cs, Co, Cr, b, c, f0):
        for m in [self]+self.children:
            m.C.Cs = Cs
            m.C.Co = Co
            m.C.Cr = Cr
            m.C.b = b
            m.C.c = c
            m.C.f0 = f0


    def collectU(self, ulist, model='do', get_mass=1):
        m_dict = {'do': self.delta_only, 'dp': self.delta_phi, 'dps': self.delta_phi_sigma}
        m = m_dict[model]
        res = []
        for U in ulist:
            xs_d = self.getDeltaXs(U)
            xo = self.delta_phi.C.X_o[10]
            xr = self.delta_phi.C.X_r[10]
            self.setDeltaConst(np.array([xs_d for i in range(4)]),
                     np.array([xo for i in range(4)]),
                     np.array([xr for i in range(4)]),
                     'xs=%.2f U = %.2f'%(xs_d, U))
            # print(wr.delta_phi.foldername)

            # m = wr.delta_only
            m.loadEos()
            if get_mass:
                mass = np.loadtxt(join(m.foldername, m.filenames['mass_crust']+'_linear'), skiprows=1)
                mmax = max(mass[:, 1])
                iM = interp1d(mass[:, 0], -mass[:, 1], kind='cubic')
                nmmax = mass[np.argmax(mass[:, 1]), 0]
                nmmax = minimize(iM, nmmax).x
            else:
                mmax = 0.
                nmmax = 0.
            lept = m.lepton_concentrations(ret_du=1)
            conc = m.concentrations()
            # plt.plot(m.nrange/m.n0, lept[2])
            # plt.plot(m.nrange/m.n0, conc[:, 1])
            # plt.show()
            # plt.plot(m.nrange/m.n0, np.abs(lept[2] - conc[:, 1]))
            # plt.show()
            # print(np.argmin(np.abs(lept[2] - conc[:, 1])))

            # n_du = m.nrange[np.argmin(np.abs(lept[2] - conc[:, 1]))]/m.n0
            if get_mass:
                n_du, m_du = m.getDuCrit()
            else:
                n_du = 0.
                m_du = 0.
            n_du /= m.n0
            _set = [0 for i in range(m.n_baryon - 2)]
            nc = [max(m.nrange/m.n0) for s in _set]
            for i, r in enumerate(m.concentrations()):
                # print(r)
                for j, s in enumerate(_set):
                    if not s:
                        if r[-len(_set)+j] > 1e-7:
                            nc[j] = m.nrange[i] / m.n0
                            _set[j] = 1

            res.append([U] + [mmax, n_du] + nc + [m_du] + [nmmax])
        np.savetxt(join(self.foldername, m.filenames['var_u']), res, fmt='%.6f')
        return np.array(res)

    def collectMbU(self, ulist, model='do'):
        res = []
        m_dict = {'do': self.delta_only, 'dp': self.delta_phi, 'dps': self.delta_phi_sigma}
        m = m_dict[model]
        for U in ulist:
            self.setDeltaPotential(U)
            mass = np.loadtxt(join(m.foldername, m.filenames['mass_crust']+'_linear'), skiprows=1)
            mmax = max(mass[:, 1])
            i_stop = np.argwhere(np.diff(mass[:, 1]) < 0)[0]
            iMb_M = interp1d(mass[:i_stop, 1], mass[:i_stop, 4], kind='cubic')
            # plt.plot(mass[:i_stop, 1], iMb_M(mass[:i_stop, 1]))
            # plt.show()
            # break
            out = iMb_M(1.249)
            res.append(out)
            print(U, iMb_M(1.249))
        np.savetxt(join(self.foldername, 'mb_u_'+model+'.dat'),
            np.array([ulist, res]).transpose())
        return np.array(ulist), np.array(res)

    def getSecondBranch(self, f_init):
        self.loadEos()
        i_start = len(self.nrange)-1
        i_end = 0


    def collectSym(self, ulist):
        res = []
        for U in ulist:
            xs_d = self.getDeltaXs(U)
            self.setDeltaConst(np.array([xs_d for i in range(4)]),
                             np.array([1. for i in range(4)]),
                             np.array([1., 1., 1., 1.]),
                             'xs=%.2f U = %.2f' % (xs_d, U))
            # print(wr.delta_sym.foldername)
            m = self.delta_sym
            m.loadEos()
            # mass = np.loadtxt(join(m.foldername, m.filenames['mass_crust']+'_linear'), skiprows=1)
            # mmax = max(mass[:, 1])
            # lept = m.lepton_concentrations(ret_du=1)
            # plt.plot(m.nrange/m.n0, m.rho)
            conc = m.concentrations()
            # plt.plot(m.nrange/m.n0, lept[2])
            # plt.plot(m.nrange/m.n0, conc[:, 1])
            # plt.show()
            # plt.plot(m.nrange/m.n0, conc)
            # plt.show()
            nc = max(m.nrange)/wr.n0
            _set = 0
            nm = max(m.nrange)/wr.n0
            for i, r in enumerate(m.concentrations()):
                if not _set:
                    if r[10] > 1e-7:
                        nc = m.nrange[i] / m.n0
                        _set = 1

            for i, r in enumerate(m.rho):
                # print(r)
                if not (r[0] < 1.):
                    nm = m.nrange[i]/m.n0
                    break
            fjump = 1.
            njump = 8.
            f = m.rho[:, 0]
            if any(np.diff(f) > 0.02):
                fjump = f[np.argmax(np.diff(f)) + 1]
                njump = m.nrange[np.argmax(np.diff(f)) + 1]/m.n0
            res.append([U, nc, nm, fjump, njump])
        np.savetxt(join(self.foldername, 'nc_sym_U.dat'), res, fmt='%.6f')


    def collectXo(self, xo_list, U, model='do', get_mass=1):
        m_dict = {'do': self.delta_only, 'dp': self.delta_phi, 'dps': self.delta_phi_sigma}
        m = m_dict[model]
        res = []
        for xo in xo_list:
            # xs_d = self.getDeltaXs(U)
            # xo = self.delta_phi.C.X_o[10]
            # xr = self.delta_phi.C.X_r[10]
            self.setDeltaXo(xo)
            self.setDeltaPotential(U)
            # self.setDeltaConst(np.array([xs_d for i in range(4)]),
            #          np.array([xo for i in range(4)]),
            #          np.array([xr for i in range(4)]),
            #          'xs=%.2f U = %.2f'%(xs_d, U))
            # print(wr.delta_phi.foldername)

            # m = wr.delta_only
            m.loadEos()
            if get_mass:
                mass = np.loadtxt(join(m.foldername, m.filenames['mass_crust']+'_linear'), skiprows=1)
                mmax = max(mass[:, 1])
                iM = interp1d(mass[:, 0], -mass[:, 1], kind='cubic')
                nmmax = mass[np.argmax(mass[:, 1]), 0]
                nmmax = minimize(iM, nmmax).x
            else:
                mmax = 0.
                nmmax = 0.
            lept = m.lepton_concentrations(ret_du=1)
            conc = m.concentrations()
            # plt.plot(m.nrange/m.n0, lept[2])
            # plt.plot(m.nrange/m.n0, conc[:, 1])
            # plt.show()
            # plt.plot(m.nrange/m.n0, np.abs(lept[2] - conc[:, 1]))
            # plt.show()
            # print(np.argmin(np.abs(lept[2] - conc[:, 1])))

            # n_du = m.nrange[np.argmin(np.abs(lept[2] - conc[:, 1]))]/m.n0
            if get_mass:
                n_du, m_du = m.getDuCrit()
            else:
                n_du = 0.
                m_du = 0.
            n_du /= m.n0
            _set = [0 for i in range(m.n_baryon - 2)]
            nc = [max(m.nrange/m.n0) for s in _set]
            for i, r in enumerate(m.concentrations()):
                # print(r)
                for j, s in enumerate(_set):
                    if not s:
                        if r[-len(_set)+j] > 1e-7:
                            nc[j] = m.nrange[i] / m.n0
                            _set[j] = 1

            res.append([xo] + [mmax, n_du] + nc + [m_du] + [nmmax])
            np.savetxt(join(self.foldername, m.filenames['var_xo']), res, fmt='%.6f')
        return np.array(res)

    def collectXr(self, xr_list, U, model='do', get_mass=1):
        m_dict = {'do': self.delta_only, 'dp': self.delta_phi, 'dps': self.delta_phi_sigma}
        m = m_dict[model]
        res = []
        for xr in xr_list:
            # xs_d = self.getDeltaXs(U)
            # xo = self.delta_phi.C.X_o[10]
            # xr = self.delta_phi.C.X_r[10]
            self.setDeltaPotential(U)
            self.setDeltaXr(xr)
            # self.setDeltaConst(np.array([xs_d for i in range(4)]),
            #          np.array([xo for i in range(4)]),
            #          np.array([xr for i in range(4)]),
            #          'xs=%.2f U = %.2f'%(xs_d, U))
            # print(wr.delta_phi.foldername)

            # m = wr.delta_only
            m.loadEos()
            if get_mass:
                mass = np.loadtxt(join(m.foldername, m.filenames['mass_crust']+'_linear'), skiprows=1)
                mmax = max(mass[:, 1])
                iM = interp1d(mass[:, 0], -mass[:, 1], kind='cubic')
                nmmax = mass[np.argmax(mass[:, 1]), 0]
                nmmax = minimize(iM, nmmax).x
            else:
                mmax = 0.
                nmmax = 0.
            lept = m.lepton_concentrations(ret_du=1)
            conc = m.concentrations()
            # plt.plot(m.nrange/m.n0, lept[2])
            # plt.plot(m.nrange/m.n0, conc[:, 1])
            # plt.show()
            # plt.plot(m.nrange/m.n0, np.abs(lept[2] - conc[:, 1]))
            # plt.show()
            # print(np.argmin(np.abs(lept[2] - conc[:, 1])))

            # n_du = m.nrange[np.argmin(np.abs(lept[2] - conc[:, 1]))]/m.n0
            if get_mass:
                n_du, m_du = m.getDuCrit()
            else:
                n_du = 0.
                m_du = 0.
            n_du /= m.n0
            _set = [0 for i in range(m.n_baryon - 2)]
            nc = [max(m.nrange/m.n0) for s in _set]
            for i, r in enumerate(m.concentrations()):
                # print(r)
                for j, s in enumerate(_set):
                    if not s:
                        if r[-len(_set)+j] > 1e-7:
                            nc[j] = m.nrange[i] / m.n0
                            _set[j] = 1

            res.append([xr] + [mmax, n_du] + nc + [m_du] + [nmmax])
            np.savetxt(join(self.foldername, m.filenames['var_xo']), res, fmt='%.6f')
        return np.array(res)


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
        self.dumpScalings()
        self.dumpBaryonParams()
        self.dumpProps()
        self.dumpEos()
        self.nucl.dumpMeff()
        self.nucl.dumpMassesCrust()
        self.nucl.dumpEos()
        for s in [self.nucl, self.sym, self.neutr]:
            s.dumpVs()
            s.dumpParts()

        mods = []
        if hyper:
            # mods = [self.hyper, self.hyper_phi, self.hyper_phi_sigma]
            mods = [self.hyper_phi, self.hyper_phi_sigma]
            # mods = [self.hyper_phi_sigma]
        for s in mods:
            s.dumpEos()
            s.dumpVs()
            try:
                s.dumpMassesCrust(npoints=100)
            except:
                pass
            s.dumpMeff()
            s.dumpEtap()

    def dumpBaryonParams(self):
        with open(join(self.foldername, 'b_params.dat'), 'w') as f:
            for i in range(8):
                f.write("X_s[%i]=%.6f"%(i, self.C.X_s[i]))
            f.write("\n")

            for i in range(8):
                f.write("X_o[%i]=%.6f"%(i, self.C.X_o[i]))
            f.write("\n")

            for i in range(8):
                f.write("X_r[%i]=%.6f"%(i, self.C.X_r[i]))
            f.write("\n")

            for i in range(8):
                f.write("X_p[%i]=%.6f"%(i, self.C.X_p[i]))
            f.write("\n")

            for i in range(8, 12):
                f.write("X_s[%i]=%.6f"%(i, self.delta_phi.C.X_s[i]))
            f.write("\n")

            for i in range(8, 12):
                f.write("X_s[%i]=%.6f"%(i, self.delta_phi.C.X_s[i]))
            f.write("\n")

            for i in range(8, 12):
                f.write("X_o[%i]=%.6f"%(i, self.delta_phi.C.X_o[i]))
            f.write("\n")

            for i in range(8, 12):
                f.write("X_r[%i]=%.6f"%(i, self.delta_phi.C.X_r[i]))
            f.write("\n")

            for i in range(8, 12):
                f.write("X_p[%i]=%.6f"%(i, self.delta_phi.C.X_p[i]))
            f.write("\n")

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
        # print(table)
        # print(join(self.foldername, 'pot_sym.dat'))
        with open(join(self.foldername, 'pot_sym.dat'), 'w') as f:
            f.write(table)
        return [arr(slist), arr(vlist)]

    def getPots(self, n):
        #We have to approach the density n from below gradually:
        steps = 10
        nrange = np.linspace(0, n, steps)
        f = [0.]
        for _n in nrange:
            n_in = np.array([_n/2, _n/2])
            f = eos.f_eq(n_in, np.array(f), 1, self.C)
        
        pots = eos.potentials(np.insert(n_in, 0, f), 5, self.C)
        V = self.m_pi*self.C.X_o[0] * pots[2]
        S = self.m_pi*(-self.C.M[0]) * (1 - self.C.phi_n(0, self.C.X_s[0]*self.C.M[0]/self.C.M[0]*f[0]))

        return [S, V]

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

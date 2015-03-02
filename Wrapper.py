import eosWrap as eos
import numpy as np
import multiprocessing as mp
import multi_progress as prog
from numpy import dtype, argmin
from scipy import interpolate, optimize
import matplotlib.pyplot as plt
import sys
import os
from tabulate import tabulate
from scipy.misc.common import derivative
from zipfile import ZipFile
from glob import glob
from os.path import join
from pylab import sqrt, pi
from AuxIntegrals import I1num, I2num, I3num, I4num
from scipy.interpolate import interp1d, interp2d
from types import *
from scipy import integrate
from scipy.optimize import root
from scipy.optimize import bisect
import inspect
from multiprocessing import Queue, Process
import time


class Wrapper():
    def __init__(self, C):
        self.C = C
        self.m_pi = 135.0
        self.m_e = 0.5/self.m_pi
        self.m_mu = 105./self.m_pi
        self.n0 = C.n0
        self.set = False
        self.driverSet = False
        self.const = 197.33*(0.16/self.n0)**(4.0/3.0)*self.m_pi**(-4.0)
        
    def stepE(self, n, last, f, lenLast, iter, C, que):
        rho = eos.stepE(n, last, f, len(last), iter, C) 
        que.put(rho)
        
    def reset(self, hyper=False, nmin=0.1, nmax = 4.0, npoints = 400,
               iter=30, verbose=0, timeout=None):
        if hyper:
            init = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        else:
            init = np.array([0.0], dtype=dtype('float'))

        self.n = np.linspace(nmin, nmax, npoints, retstep=1, endpoint=0)
        print self.n[1]/self.n0
        self.n = self.n[0]
        last = init
        self.rho=[]


#         p = mp.Pool()
        writer = prog.Writer((0,0))
        pb = prog.ProgressBar(fd=writer)
        pb.start()
        self.mu_e = []

        sp = 1 + self.C.sprime
        f = np.array([0.0 for i in range(sp)])
        for _i, i in enumerate(self.n):
            if abs(i - 3.53) < 1e-2:
                pass

            pb.update(100 * _i / len(self.n))
            if timeout is None:
                rho = eos.stepE(i, last, f, len(last), iter, self.C)
            else:
                queue = Queue()
                p = Process(target=self.stepE, args=(i, last, f, len(last),
                                                  iter, self.C, queue))
                p.start()
                p.join(timeout)
                if p.is_alive():
                    p.terminate()
                    print "timeout reached"
                    self.rho = np.ascontiguousarray(self.rho[:])
                    self.n = self.n[: self.rho.shape[0]]
                    _E = map(lambda z: eos.E(z, self.C), self.rho)
                    self.E = np.ascontiguousarray(np.array(_E[:]))
                    self.P = np.ascontiguousarray(self.P_chem(self.rho))
                    self.set=1
                    return
#                     p.join() #Needed???
                rho = queue.get(timeout=1)
                print rho
#                 exit()
    
            self.rho.append(rho.copy())
            self.rho[_i]=np.insert(self.rho[_i],0, i - np.sum(rho))
            f = eos.f_eq(self.rho[_i], f, sp, self.C)
            if verbose:
                print f, rho
                print self.C.eta_r(f[0]), self.C.Cr/self.C.eta_r(f[0])
            for j in range(sp):
                self.rho[_i]=np.insert(self.rho[_i], j, f[j])

            self.mu_e.append(eos.mu(self.rho[_i], 1, self.C) - eos.mu(self.rho[_i], 2, self.C))
            last = np.array(rho)
#             temp = self.rho[_i]
#             temp[3] = 0.0
#             print eos.mu(self.rho[_i], 1, self.C), eos.mu(temp, 3, self.C)

        self.rho=np.array(self.rho)
        _E = map(lambda z: eos.E(z, self.C), self.rho)
        dn = self.n[1] - self.n[0]
#         self.P = self.n[1:]*np.diff(_E)/dn - _E[1:]
#         self.P = np.insert(self.P, 0, 0.)
#         self.P = self.n*self.diffArray(_E, dn) - _E[:]
#         self.P = np.insert(self.P, 0, 0.)
        self.E = np.ascontiguousarray(np.array(_E[:]))

#         self.P = self.n*self.diffArray(self.E, dn) - self.E[:]

        self.n = np.ascontiguousarray(self.n[:])
        self.P = np.ascontiguousarray(self.P_chem(self.rho))
        self.rho = np.ascontiguousarray(self.rho[:])
        self.set = True

    def EPN(self):
        if self.set:
            return self.E*self.m_pi**4, self.P*self.m_pi**4, self.n
        else:
            print 'Wrapper is not set!'

    def setDriver(self):
        if not self.set:
            print 'Wrapper is not set!'
            return

        self.dr = eos.KVDriver()
        self.dr.set(self.E*self.m_pi**4, self.P*self.m_pi**4, self.n)
        self.driverSet = True

    def solve(self, f0=0.195, E0=-16.0, K0=275.0, J0=32.0,
              iter=1000, mu_scale = 1):
        sprime_old = self.C.sprime
        self.C.sprime = 0
        res = eos.solve(self.n0, E0, f0, K0, J0, self.C, iter, mu_scale)
        print res
        if mu_scale < 50000 and res == 5:
            self.solve(f0, E0, K0, J0, iter, 5*mu_scale)
        self.C.sprime = sprime_old

    def stars(self, nmin = 0.5, nmax = 4.0, npoints=20):
        if not self.set:
            print 'Wrapper is not set!'
            return

        if not self.driverSet:
            print 'Driver is not set!'
            return

        dimRes = 3

        self.dr = eos.KVDriver()
        self.dr.set(self.E, self.P, self.n)

        self.n_star = np.linspace(nmin, nmax, npoints)
        res = np.array(map(lambda z: eos.star_crust2(z, dimRes, self.dr,
                                                      0.8*self.n0),
                            self.n_star))
        if dimRes == 3:
            return (self.n_star, res[:,0], res[:,1],
                    self.m_pi**4 * self.C.M[0] * res[:,2])


    def Esymm(self, n, f=0.0, ret_f=False):
        res = []
        flist = []
        for z in n:
#             print z
            f, = eos.f_eq(np.array([z/2,z/2]), np.array([f]), 1, self.C)
#             res.append((eos.E(np.array([f, z/2, z/2]), self.C))/z - self.C.M[0])
            res.append((eos.E(np.array([f, z/2, z/2]), self.C)))
            flist.append(f)
        if not ret_f:
            return np.array(res)
        else:
            return np.array(res), np.array(flist)

    def Eneutr(self, n, f=0.0, ret_f=False):
        res = []
        flist = []
        for z in n:
            f, = eos.f_eq(np.array([z,0.]), np.array([f]), 1, self.C)
#             res.append((eos.E(np.array([f, z/2, z/2]), self.C))/z - self.C.M[0])
            res.append((eos._E(np.array([f, z, 0.]), self.C)))
            flist.append(f)
        if not ret_f:
            return np.array(res)
        else:
            return np.array(res), np.array(flist)


    def Psymm(self, n):
        nlist = []
        f = 0.
        for _n in n:
            f,= eos.f_eq(np.array([_n/2, _n/2]), np.array([f]), 1, self.C)
            nlist.append(np.array([f, _n/2, _n/2]))


#         res = n[1:]*np.diff(E)/(n[1]-n[0]) - E[1:]
#         res = np.insert(res, 0, 0.)
#         res = n[:]*self.diffArray(E,n[1]-n[0]) - E[:]
        res = self.P_chem(nlist)
#         res = np.insert(res, 0, 0.)
        return self.m_pi**4 * self.const * res

    def P_N(self, n):
        nlist = []
        f = 0.
        for _n in n:
            f,= eos.f_eq(np.array([_n, 0.]), np.array([f]), 1, self.C)
            nlist.append(np.array([f, _n, 0.]))


#         res = n[1:]*np.diff(E)/(n[1]-n[0]) - E[1:]
#         res = np.insert(res, 0, 0.)
#         res = n[:]*self.diffArray(E,n[1]-n[0]) - E[:]
        res = self.P_chem(nlist, leptons=False)
#         res = np.insert(res, 0, 0.)
        return self.m_pi**4 * self.const * res

    def star_print(self, z):
        return eos.star_crust(z, 3, self.dr, 1e-11)

    def stars_crust(self, ncut_crust=0.6, ncut_eos = 0.9, inter='linear', nmin = .4, nmax = 4.0,
                     npoints=50, crust="crust.dat", show=False, crustName=None):
        neos= 1000
        if self.C.Hyper:
            neos = 1000
        self.reset(nmin = 0., hyper=self.C.Hyper, nmax=nmax, npoints=neos, timeout=5)
        E, P, N = self.EPN()
        N = N/self.n0
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

        plist = np.append(p[:],P[i_n_eos:(i_n_eos+n_eos)])
        elist = np.append(e[:],E[i_n_eos:(i_n_eos+n_eos)])
        nlist = np.append(n[:],N[i_n_eos:(i_n_eos+n_eos)])


        print i_n_eos
        print P[i_n_eos:(i_n_eos+n_eos)]
        print N[i_n_eos:(i_n_eos+n_eos)]
#         exit()
        iP = interpolate.interp1d(nlist, plist, kind=inter)
        iE = interpolate.interp1d(nlist, elist, kind=inter)

        gamma = 1./4.
        iN = np.linspace(1e-10**gamma, ncut_eos**gamma, 10000)
        iN = iN**(1./gamma)

        crust_p = np.array(map(iP, iN))
        crust_e = np.array(map(iE, iN))

#         finalE = np.append(crust_e, E[i_n_eos+n_eos:])
#         finalP = np.append(crust_p, P[i_n_eos + n_eos:])
#         finalN = np.append(iN, N[i_n_eos+n_eos:])

        finalE = np.append(crust_e, E[i_n_eos+n_eos:])/self.m_pi**4
        finalP = np.append(crust_p, P[i_n_eos + n_eos:])/self.m_pi**4
        finalN = np.append(iN, N[i_n_eos+n_eos:])

        if show:
            lines = plt.plot(finalN, finalP*self.m_pi**4, n, p, N, P)
            plt.xlabel(r'$n/n_0$', fontsize=18)
            plt.xlim([0, 2])
            plt.ylabel(r'$P \, [MeV^4] $', fontsize=18)
            plt.legend(lines, ['interpolated', 'crust', 'RMF'], loc=0)
            plt.show()

        if crustName is not None:
            tab = np.array([finalN, finalE, finalP]).transpose()
            table = tabulate(tab, ['n/n0', 'E','P'], tablefmt='plain')
            with open(crustName, 'w') as f:
                f.write(table)

        self.dr = eos.KVDriver()
        self.dr.set(finalE, finalP, finalN*self.n0)
        nstar = np.linspace(nmin, nmax, npoints)

        MR = []
        Mgrav = []
        for _n in nstar:
            MR.append(eos.star_crust2(_n, 3, self.dr, 1e-11))
            lastN = self.dr.getLastN(self.dr.nSize)[:-1]
            lastR = self.dr.getLastR(self.dr.nSize)[:-1]
            lastM = self.dr.getLastM(self.dr.nSize)[:-1]
            dx = lastR[1] - lastR[0]

            grav_mult = []
            for i, r in enumerate(lastR):
                grav_mult.append( r**2 / sqrt(1 - 2 * 1.4677 * lastM[i] / r))

            grav_mult = np.array(grav_mult)
            res = np.multiply(lastN, grav_mult)
            Mgrav.append((0.0004898007281478712)*np.trapz(res, dx=dx))

        MR = np.array(MR)
        Mgrav = np.array(Mgrav)
        print MR
        print Mgrav
#         return nstar, MR[:, 0], MR[:,1], self.m_pi**4*self.C.M[0]*MR[:,2]
        return (nstar, MR[:, 0], MR[:,1], 931.5/self.m_pi*MR[:,2],
                931.5/self.m_pi * Mgrav)

    def stars_crust_hyper(self, ncut_crust=0.6, ncut_eos = 0.8,
                         inter='linear', nmin = .4, nmax = 4.0,
                         npoints=100, crust="crust.dat", show=False):
        hyp_old = self.C.Hyper
        self.C.Hyper=1
        res = self.stars_crust(ncut_crust=ncut_crust, ncut_eos = ncut_eos,
                         inter=inter, nmin = nmin, nmax = nmax,
                         npoints=npoints, crust=crust, show=show)
        self.C.Hyper=hyp_old
        return res

    def concentrations(self):
        if not self.set:
            print 'Concentrations: wrapper is not set! Pass.'
            return

        rho = []
        sp = 1 + self.C.sprime
        for _r in self.rho:
            if np.sum(_r[sp:]) > 0:
                rho.append(_r[sp:]/np.sum(_r[sp:]))
            else:
                rho.append([1.]+[0. for i in _r[sp+1:]])

        return np.array(rho)

    def dumpEos(self, folderName):
        if not os.path.exists(folderName):
            os.makedirs(folderName)

#         if not self.set:
        self.reset(nmin = 0., nmax = 8. * self.n0, npoints=800)

        Esymm = []
        f = 1e-6
        flist_symm = []
        nlist_symm = []
        for n in self.n:
#             print f
#             print eos.EBind(np.array([f, n/2, n/2]), self.C)
            f, = eos.f_eq(np.array([n/2, n/2]), np.array([f]), 1, self.C)
            flist_symm.append(f)
            nlist_symm.append(np.array([n/2, n/2]))
            Esymm.append(eos.EBind(np.array([f, n/2, n/2]), self.C))

        Esymm[0] = 0.
#         n_symm = np.array(nlist_symm)
        Psymm = self.Psymm(self.n)

        EN = []
        eps_N=[]
        f = 1e-6
        flist_n = []

        for n in self.n:
#             print f
#             print eos.EBind(np.array([f, n/2, n/2]), self.C)
            f, = eos.f_eq(np.array([n, 0.]), np.array([f]), 1, self.C)
            flist_n.append(f)
            EN.append(eos.EBind(np.array([f, n, 0.]), self.C))
            eps_N.append(eos._E(np.array([f, n, 0]), self.C))

        EN[0] = 0.

#         PN = self.n[1:]*np.diff(EN)/(self.n[1]-self.n[0]) - EN[1:]
#         PN = np.insert(PN, 0, 0.)
        eps_N = np.array(eps_N)
        PN = self.P_N(self.n)
#         PN = self.n[:]*self.diffArray(eps_N, self.n[1]-self.n[0]) - eps_N[:]
#
#         PN *= self.m_pi**4 * self.const

        rho = self.concentrations()

        table = []
        
        x_du=[]
        l_n_e = []
        l_n_mu = []
        for _mu in self.mu_e:
            n_e = 0
            n_mu = 0
            if _mu > self.m_e:
                n_e += (_mu**2 - self.m_e**2)**(1.5)/(3*pi**2)
            if _mu > self.m_mu:
                n_mu += (_mu**2 - self.m_mu**2)**(1.5)/(3*pi**2)
            
            l_n_e.append(n_e)
            l_n_mu.append(n_mu)
            
            if n_e + n_mu > 0:
                x_du.append(1./(1 + (1 + (n_e/(n_e + n_mu))**(1./3))**3))
            else:
                x_du.append(0.11)
        
        for i, _n in enumerate(self.n):
            table.append([_n/self.n0, Esymm[i], Psymm[i], flist_symm[i],
                         self.E[i], self.P[i], rho[i,1], self.rho[i, 0],
                         EN[i], PN[i], flist_n[i], self.mu_e[i]*self.m_pi,
                         x_du[i]])

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

        with open(os.path.join(folderName, 'eos.dat'), 'w') as f:
            f.write(tab)

        return

    def dumpMasses(self, folderName=None):
        if folderName is not None:
            if not os.path.exists(folderName):
                os.makedirs(folderName)

        if not self.set:
            self.reset(nmin = 0., nmax = 8. * self.n0, npoints=800)

        if not self.driverSet:
            self.setDriver()

        n, m, r, mb = self.stars(npoints = 100)
        if folderName is not None:
            table = np.array([n/self.n0, m, r, mb]).transpose()
            f = open(os.path.join(folderName, 'masses.dat'), 'w')
            tab = tabulate(table, ['n/n_0',
                                   'M [M_{sun}]',
                                   'R [km]',
                                   'M_B [M_sun]'],
                           tablefmt='plain')
            f.write(tab)
            f.close()
        return n, m, r

    def dumpHyper(self, folderName, npoints=400, verbose=False, filename=''):
        if not os.path.exists(folderName):
            os.makedirs(folderName)
    
        #SCALINGS(f) no phi 
        frange = np.linspace(0, 1, 100)
        phi_sc = map(self.C.eta_p, frange)
        etap_tab = np.array([frange, phi_sc]).transpose()
        etap_table = tabulate(etap_tab, ['f', 'eta_p'], tablefmt='plain')
        with open(join(folderName, 'etap.dat'), 'w') as f:
            f.write(etap_table)
        
        tab = []
        for _f in frange:
            _tab = []
            _tab.append(_f)
            for i in xrange(2, 8):
                _tab.append(self.C.Xs(i, _f)/self.C.X_s[i])
            tab.append(_tab)
        etap_tab = np.array(tab)#.transpose()
        etap_table = tabulate(etap_tab, ['f', 'L','S-','S0','S+','X-','X0'], tablefmt='plain')
        with open(join(folderName, 'xi_s.dat'), 'w') as f:
            f.write(etap_table)
            
        #RESET no phi
        self.C.phi_meson = 0
        self.reset(hyper = 1, nmin = 0., nmax = 8. * self.n0,
                    npoints=npoints, iter=100, verbose=verbose, timeout=6)
        
        #MASSES no phi
        self.setDriver()
        n, m, r, mg1, mg2 = self.stars_crust_hyper(npoints=100)
        
        mtab = np.array([n/self.n0, m, r]).transpose()

        mtable = tabulate(mtab, ['n/n_0','M/M_sun', 'R [km]'], tablefmt='plain')

        with open(join(folderName, 'mass_hyper'+filename+'.dat'), 'w') as f:
            f.write(mtable)
            
        #CONCENTRATIONS no phi
        rho = self.concentrations()
        print rho
        table = []
        for i, _n in enumerate(self.n):
            table.append(np.concatenate(([_n/self.n0, self.E[i],
                                           self.P[i]], rho[i], [self.rho[i,0]] )))

        tab = tabulate(table, ['n/n_0',
                               'E_NS [m_\pi^4]',
                               'P_NS, [m_\pi^4]',
                               'x_n',
                               'x_p',
                               'x_Lambda',
                               'x_Sigma-', 'x_Sigma0', 'x_Sigma+',
                               'x_Xi-', 'x_Xi0', 'f(n)'],
                       tablefmt='plain')

        with open(os.path.join(folderName, 'hyper'+filename+'.dat'), 'w') as f:
            f.write(tab)

        #EFF. MASSES no phi

        tab = []
        for i, n in enumerate(self.n):
            f = self.rho[i, 0]
            tab.append([n/self.n0,
                        self.C.phi_n(0, f),
                        self.C.phi_n(2, self.C.Xs(2,f) * f),
                        self.C.phi_n(3, self.C.Xs(3,f) * f),
                        self.C.phi_n(6, self.C.Xs(6,f) * f)])

        table = tabulate(tab, ['n/n0', 'N', 'Lambda', 'Sigma', 'Xi'],
                         tablefmt='plain')

        with open(join(folderName, 'hyper_meff'+filename+'.dat'), 'w') as f:
            f.write(table)

        #SOUND SPEED no phi
        
        vsNs = np.gradient(self.P)/np.gradient(self.E)
        tab = np.array([self.n[:]/self.n0, vsNs]).transpose()
        table = tabulate(tab, ['n/n0', 'vs^2'], tablefmt='plain')
        
        if folderName is not None:
            with open(join(folderName, 'vsHyper.dat'), 'w') as f:
                f.write(table)            
        
        #SCALINGS 2 no phi
        frange = np.linspace(0, 1, 100)
        phi_sc = map(self.C.eta_p, frange)
        chi_p = sqrt((1 - frange)**2/np.array(phi_sc))
        etap_tab = np.array([frange, phi_sc]).transpose()
        etap_table = tabulate(etap_tab, ['f', 'eta_p'], tablefmt='plain')
        with open(join(folderName, 'etap.dat'), 'w') as f:
            f.write(etap_table)
#         print chi_p
        chi_tab = np.array([frange, chi_p]).transpose()
#         print chi_tab
        chi_table = tabulate(chi_tab, ['f', 'chi_p'], tablefmt='plain')
        with open(join(folderName, 'chi_p.dat'), 'w') as f:
            f.write(chi_table)
        
        tab = []
        for _f in frange:
            _tab = []
            _tab.append(_f)
            for i in xrange(2, 8):
                _tab.append(self.C.Xs(i, _f)/self.C.X_s[i])
            tab.append(_tab)
        etap_tab = np.array(tab)#.transpose()
        etap_table = tabulate(etap_tab, ['f', 'L','S-','S0','S+','X-','X0'], tablefmt='plain')
        with open(join(folderName, 'xi_s.dat'), 'w') as f:
            f.write(etap_table)
            
        frange = self.rho[:, 0]
        phi_sc = map(self.C.eta_p, frange)
        etap_tab = np.array([self.n/self.n0, phi_sc]).transpose()
        etap_table = tabulate(etap_tab, ['n/n0', 'eta_p'], tablefmt='plain')
        with open(join(folderName, 'etap_N.dat'), 'w') as f:
            f.write(etap_table)
        
        chi_p = sqrt((1 - frange)**2/np.array(phi_sc))
        chi_tab = np.array([self.n/self.n0, chi_p]).transpose()
        chi_table = tabulate(chi_tab, ['n/n0', 'chi_p'], tablefmt='plain')
        with open(join(folderName, 'chi_p_N.dat'), 'w') as f:
            f.write(chi_table)
        
        tab = []
        for j, _f in enumerate(frange):
            _tab = []
            _tab.append(self.n[j]/self.n0)
            for i in xrange(2, 8):
                _tab.append(self.C.Xs(i, _f)/self.C.X_s[i])
            tab.append(_tab)
        etap_tab = np.array(tab)#.transpose()
        etap_table = tabulate(etap_tab, ['n/n0', 'L','S-','S0','S+','X-','X0'], tablefmt='plain')
        with open(join(folderName, 'xi_s_N.dat'), 'w') as f:
            f.write(etap_table)
        
        #RESET phi
        self.C.phi_meson = 1
        self.reset(hyper = 1, nmin = 0., nmax = 8. * self.n0,
                    npoints=npoints, iter=100)
 
        #MASSES phi
        self.setDriver()
        n, m, r, mg1, mg2 = self.stars_crust_hyper(npoints=100)

        mtab = np.array([n/self.n0, m, r]).transpose()

        mtable = tabulate(mtab, ['n/n_0','M/M_sun', 'R [km]'], tablefmt='plain')

        with open(join(folderName, 'mass_hyper_phi'+filename+'.dat'), 'w') as f:
            f.write(mtable)

        #CONCENTRATIONS phi
        rho = self.concentrations()
        table = []
        for i, _n in enumerate(self.n):
            table.append(np.concatenate(([_n/self.n0, self.E[i],
                                           self.P[i]], rho[i], [self.rho[i, 0]])))

        tab = tabulate(table, ['n/n_0',
                               'E_NS [m_\pi^4]',
                               'P_NS, [m_\pi^4]',
                               'x_n',
                               'x_p',
                               'x_Lambda',
                               'x_Sigma-', 'x_Sigma0', 'x_Sigma+',
                               'x_Xi-', 'x_Xi0', 'f'],
                       tablefmt='plain')

        with open(os.path.join(folderName, 'hyper_phi'+filename+'.dat'), 'w') as f:
            f.write(tab)

        #EFF. MASSES phi
        tab = []
        for i, n in enumerate(self.n):
            f = self.rho[i, 0]
            tab.append([n/self.n0,
                        self.C.phi_n(0, f),
                        self.C.phi_n(2, self.C.Xs(2,f) * f),
                        self.C.phi_n(3, self.C.Xs(3,f) * f),
                        self.C.phi_n(6, self.C.Xs(6,f) * f)])

        table = tabulate(tab, ['n/n0', 'N', 'Lambda', 'Sigma', 'Xi'],
                         tablefmt='plain')

        with open(join(folderName, 'hyper_meff_phi'+filename+'.dat'), 'w') as f:
            f.write(table)

        #SOUND SPEED phi
        vsNs = np.gradient(self.P)/np.gradient(self.E)
        tab = np.array([self.n[:]/self.n0, vsNs]).transpose()
        table = tabulate(tab, ['n/n0', 'vs^2'], tablefmt='plain')
        
        if folderName is not None:
            with open(join(folderName, 'vsHyperPhi.dat'), 'w') as f:
                f.write(table)

        #SCALINGS 2 phi
        frange = self.rho[:, 0]
        phi_sc = map(self.C.eta_p, frange)
        etap_tab = np.array([self.n/self.n0, phi_sc]).transpose()
        etap_table = tabulate(etap_tab, ['n/n0', 'eta_p'], tablefmt='plain')
        with open(join(folderName, 'etap_N_phi.dat'), 'w') as f:
            f.write(etap_table)
        
        chi_p = sqrt((1 - frange)**2/np.array(phi_sc))
        chi_tab = np.array([self.n/self.n0, chi_p]).transpose()
        chi_table = tabulate(chi_tab, ['n/n0', 'chi_p'], tablefmt='plain')
        with open(join(folderName, 'chi_p_N.dat'), 'w') as f:
            f.write(chi_table)
        
        tab = []
        for j, _f in enumerate(frange):
            _tab = []
            _tab.append(self.n[j]/self.n0)
            for i in xrange(2, 8):
                _tab.append(self.C.Xs(i, _f)/self.C.X_s[i])
            tab.append(_tab)
        etap_tab = np.array(tab)#.transpose()
        etap_table = tabulate(etap_tab, ['n/n0', 'L','S-','S0','S+','X-','X0'], tablefmt='plain')
        with open(join(folderName, 'xi_s_N_phi.dat'), 'w') as f:
            f.write(etap_table)   

        self.C.phi_meson = 0
        

    def dumpAll(self, folderName,filename='data.zip'):
        self.dumpEos(folderName)
        n, m ,r = self.dumpMasses(folderName)
        mmax = max(m)
        rho = self.concentrations()
        print rho[:, 1]
        print abs(rho[:,1] - [0.14 for i in rho[:,1]])
        print np.argmin(abs(rho[:,1] - [0.14 for i in rho[:,1]]))
        ndu = self.n[np.argmin(abs(rho[:,1] - [0.14 for i in rho[:,1]]))]
        mdu = m[np.argmin(abs(n - [ndu for i in n]))]
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
                  ['M_max [M_sun]', mmax],
                  ['M_DU [M_sun]', mdu],
                  ['n_DU [n_0]', ndu/self.n0],
                  ['n_0 [fm^-3]', self.n0 * (135/197.33)**3]
                  ]
        Params = tabulate(tabParams,floatfmt='.10f',tablefmt='plain')
        Part = tabulate(tabPart,tablefmt='plain',floatfmt='.10f')
        Res = tabulate(tabRes,tablefmt='plain',floatfmt='.10f')
        
#         cPar = ''
#         for member in inspect.getmembers(self.C, inspect.isdatadescriptor):
#             print member
#         
#         exit()
        f = open(os.path.join(folderName, 'props.dat'), 'w')
        f.write(Params + '\n' + Part + '\n' + Res + '\n')
        f.close()
        self.dumpScalings(folderName)
#         self.dumpHyper(folderName)
        with ZipFile(join(folderName, filename), 'w') as zip:
            for f in glob(join(folderName, '*.dat')):
                zip.write(f)



    def dumpScalings(self, folderName):
        if not self.set:
            self.reset(nmin = 0., nmax = 4., npoints = 1000)

        tabNS = []
        for i, _r in enumerate(self.rho):
            f = _r[0]
            C = self.C
            tabNS.append([self.n[i]/self.n0,
                          f,
                          C.eta_s(f) + 2*C.Cs*C.U(f)/(C.M[0]**4 * f**2),
                          C.eta_o(f),
                          C.eta_r(f),
                          C.phi_n(0, f),
                          C.U(f)])

        tabF = []
        for f in np.linspace(0., 1., 100., endpoint=False):
            C = self.C
            tabF.append([f, C.eta_s(f) + 2*C.Cs*C.U(f)/(C.M[0]**4 * f**2),
                          C.eta_o(f), C.eta_r(f),
                          C.eta_p(f), C.phi_n(0, f), C.U(f)])

        tableNS = tabulate(tabNS, headers=['n/n_0',
                                           'f(n)',
                                           'eta_sigma',
                                           'eta_omega',
                                           'eta_rho',
                                           'Phi_N',
                                           'U'],tablefmt='plain', floatfmt='.10f')

        tableF = tabulate(tabF, headers=['f',
                                         'eta_sigma',
                                         'eta_omega',
                                         'eta_rho',
                                         'eta_phi',
                                         'Phi_N',
                                         'U'],tablefmt='plain', floatfmt='.10f')
        if folderName is not None:
            with open(os.path.join(folderName, 'scalingsNS.dat'), 'w') as f:
                f.write(tableNS)

            with open(os.path.join(folderName, 'scalingsF.dat'), 'w') as f:
                f.write(tableF)

        return np.array(tabNS)


    def dumpMassesCrustHyper(self, folderName ,ncut_crust=.6, ncut_eos=.8):
        if not os.path.exists(folderName):
            os.makedirs(folderName)
        
        hyper_old = self.C.Hyper
        self.C.Hyper = 1
        
        n, m, r, mb1, mb2 = self.stars_crust_hyper(ncut_crust=ncut_crust,
                                   ncut_eos=ncut_eos,npoints = 100, show=0)
        print m
        print mb1
        table = np.array([n/self.n0, m, r, mb1, mb2]).transpose()
        f = open(os.path.join(folderName, 'masses_crust_hyper.dat'), 'w')
        tab = tabulate(table, ['n/n_0','M [M_{sun}]',
                                'R [km]','M_B_rect [M_{sun}]',
                                'M_B_trap [M_{sun}]'],
                       tablefmt='plain')
        f.write(tab)
        f.close()
        self.C.Hyper = hyper_old
        
        return n, m, r, mb1, mb2
    
    def dumpMassesCrust(self, folderName ,ncut_crust=.6, ncut_eos=.8, hyper=False, show=False, inter='linear'):
        if not os.path.exists(folderName):
            os.makedirs(folderName)
        
        self.C.hyper = hyper
        
        crustName = join(folderName, 'crustEos.dat')
        n, m, r, mb1, mb2 = self.stars_crust(ncut_crust=ncut_crust,
                                   ncut_eos=ncut_eos,npoints = 100, show=show, inter=inter, crustName=crustName)
        print m
        print mb1
        table = np.array([n/self.n0, m, r, mb1, mb2]).transpose()
        f = open(os.path.join(folderName, 'masses_crust.dat'), 'w')
        tab = tabulate(table, ['n/n_0','M [M_{sun}]',
                                'R [km]','M_B_rect [M_{sun}]',
                                'M_B_trap [M_{sun}]'],
                       tablefmt='plain')
        f.write(tab)
        f.close()
        return n, m, r, mb1, mb2

    def dumpLandauParams(self, folderName):
        if not os.path.exists(folderName):
            os.makedirs(folderName)

        n_f0 = np.linspace(0., 8*self.n0, 400, endpoint=0)
        f0 = self.f0(n_f0)
        tab_f0 = np.array([n_f0/self.n0, f0]).transpose()
        table_f0 = tabulate(tab_f0, ['n/n_0', 'F_0'], tablefmt='plain')
        with open(join(folderName, 'f0.dat'), 'w') as f:
            f.write(table_f0)

        f1 = self.f1(n_f0)
        tab_f1 = np.array([n_f0/self.n0, f1]).transpose()
        table_f1 = tabulate(tab_f1, ['n/n_0', 'F_1'], tablefmt='plain')

        with open(join(folderName, 'f1.dat'), 'w') as f:
            f.write(table_f1)

        return n_f0, f0, f1
    
    def dumpLandauParamsNM(self, folderName):
        if not os.path.exists(folderName):
            os.makedirs(folderName)

        n_f0 = np.linspace(0., 8*self.n0, 400, endpoint=0)
        f0 = self.f0_nm(n_f0)
        tab_f0 = np.array([n_f0/self.n0, f0]).transpose()
        table_f0 = tabulate(tab_f0, ['n/n_0', 'F_0'], tablefmt='plain')
        with open(join(folderName, 'f0_nm.dat'), 'w') as f:
            f.write(table_f0)

        f1 = self.f1_nm(n_f0)
        tab_f1 = np.array([n_f0/self.n0, f1]).transpose()
        table_f1 = tabulate(tab_f1, ['n/n_0', 'F_1'], tablefmt='plain')

        with open(join(folderName, 'f1_nm.dat'), 'w') as f:
            f.write(table_f1)

        return n_f0, f0, f1

    def dumpPotentialsNS(self, folderName, show=False):
        if not os.path.exists(folderName):
            os.makedirs(folderName)

        self.reset()

        pots = []

        for i, _r in enumerate(self.rho):
            pots.append(self.n[i] + eos.potentials(_r, 5, self.C))

        print pots

        if show:
            plt.plot(self.n/self.n0, pots)
#             plt.plot(self.n/self.n0, np.sum(pots, axis=1)*135.)
            plt.xlabel(r'$n/n_0$')
            plt.ylabel(r'$Terms \, [m_pi^4]$')
            plt.show()

        table = tabulate(pots, ['n/n_0', 'sigma', 'sigma^*',
                                 'omega', 'rho', 'phi'],
                         tablefmt='plain')

        with open(join(folderName,'potsNS.dat'), 'w') as f:
            f.write(table)

        return pots

    def dumpHinHPotentials(self):
        sp = 1 + self.C.sprime
        pots = []
        n = np.linspace(0., 8.*self.n0, 400)
        for k in [2, 3, 7]:
            f = np.array([0.0 for i in range(sp)])
            pot = []
            for _n in n:
                n_f = np.array([0.0 for i in range(8)])
                n_f[k] = _n
                f = eos.f_eq(n_f, f, sp, self.C)
                print 'f=', f
                print 'n_f=', n_f
                n_mu = np.insert(n_f, 0, f)
                print 'insert=', n_mu
                pot.append(eos.mu(n_mu, k+sp, self.C)-self.C.M[k])
            print pot
            pots.append(pot)
        print pots
        pots=np.array(pots)
        print pots.shape, n.shape
        plt.plot(n/self.n0, pots.transpose())
        plt.show()

    def dumpMeffHyper(self, folderName):
        for phi_meson, suff in enumerate(['', '_phi']):
            self.C.phi_meson = phi_meson
            self.reset(hyper=1, nmin=0., nmax=8*self.n0, npoints=400)
            tab = []
            for i, n in enumerate(self.n):
                f = self.rho[i, 0]
                tab.append([n,
                            self.C.phi_n(0, f),
                            self.C.phi_n(2, self.C.Xs(2,f) * f),
                            self.C.phi_n(3, self.C.Xs(3,f) * f),
                            self.C.phi_n(6, self.C.Xs(6,f) * f)])

            table = tabulate(tab, ['n/n0', 'N', 'Lambda', 'Sigma', 'Xi'],
                             tablefmt='plain')

            with open(join(folderName, 'hyper_meff'+suff+'.dat')) as f:
                f.write(table)

    def testPodsiedlowski(self, n_crust=0.6, n_eos=0.9, folderName=None):
        n, m, r, mb1, mb2 = self.stars_crust(ncut_crust=n_crust, ncut_eos=n_eos)
        iM = interpolate.interp1d(n, m)
        iMb1 = interpolate.interp1d(n, mb1)
        iMb2 = interpolate.interp1d(n, mb2)

        iN = np.linspace(n[0], n[-1], 100)
        plt.plot(iMb1(iN), iM(iN))
        rect = plt.Rectangle([1.366, 1.248], 0.009, 0.002)
        plt.gca().add_patch(rect)
        plt.xlim([1.32, 1.39])
        plt.ylim([1.225, 1.275])
        plt.xlabel(r'$M_B/M_{\odot}$')
        plt.ylabel(r'$M_G/M_{\odot}$')
        if folderName is not None:
            if not os.path.exists(folderName):
                os.mkdir(folderName)
    
            plt.savefig(join(folderName,'podsiedlowski.pdf'))
        plt.show()

    def testDanielewicz(self):
        UKlahnY = []
        UKlahnX = []
        with open('/home/const/workspace2/swigEosWrapper/klahnUpper', 'r') as f:
            for line in f:
                _n, p = line.strip().split()
                UKlahnY.append(float(p))
                UKlahnX.append(float(_n)/0.16)
        plt.semilogy(UKlahnX, UKlahnY, c = 'grey')

        LKlahnX = []
        LKlahnY = []
        with open('/home/const/workspace2/swigEosWrapper/klahnLower', 'r') as f:
            for line in f:
                _n, p = line.strip().split()
                LKlahnY.append(float(p))
                LKlahnX.append(float(_n)/0.16)
        plt.semilogy(LKlahnX, LKlahnY, c = 'grey')

        n_p = np.linspace(0.0, 4.0, 800)
        plt.plot(n_p[:]/self.n0, self.Psymm(n_p))
        plt.show()

    def testDU(self, folderName=None):
        if not self.set:
            self.reset(hyper=0, npoints=1000)
        rho = self.concentrations()
        n, m, r = self.dumpMasses(folderName)
        ndu = self.n[np.argmin(abs(rho[:,1] - [0.14 for i in rho[:,1]]))]
        mdu = m[np.argmin(abs(n - [ndu for i in n]))]
        print 'n_DU = %.2f n_0, M_DU = %.2f M_sun'%(ndu/self.n0, mdu)
        return 'n_DU = %.2f n_0, M_DU = %.2f M_sun'%(ndu/self.n0, mdu)

    def testHyperBind(self, show=True):
        n = np.linspace(0., 4.*self.n0, 100, endpoint=0)
        f = 0.
        self.C.sprime = 0
        EL = []
        ES = []
        EX = []
        for _n in n:
            n_f = np.array([_n/2, _n/2]+[0.0 for i in range(6)])
            f, = eos.f_eq(n_f, np.array([f]), 1, self.C)
            print np.array([f]), n_f
            n_mu = np.insert(n_f, 0, f)
            EL.append(eos.mu(n_mu, 2+1, self.C) - self.C.M[2])
            ES.append(eos.mu(n_mu, 3+1, self.C) - self.C.M[3])
            EX.append(eos.mu(n_mu, 6+1, self.C) - self.C.M[6])

        print EL
        print ES
        print EX

        EL = self.m_pi*np.array(EL)
        ES = self.m_pi*np.array(ES)
        EX = self.m_pi*np.array(EX)


        if show:
            plt.ylim([-50, 50])
            plt.plot(n/self.n0, EL, n/self.n0, ES, n/self.n0, EX)
            plt.show()

        return EL, ES, EX


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

    def ESbind(self, n):
        return self.m_pi*(self.Esymm(n)/n - self.C.M[0])
    
    def _f0(self, n, f, multiply=True):
        C = self.C
        mn = C.M[0]
        pf = eos.p_f(n/2)
        f, = eos.f_eq(np.array([n/2, n/2]), np.array([f]), 1, C)

        meff = mn * C.phi_n(0, f)
        res = C.Co/(mn**2 * C.eta_o(f))
        dPhi = derivative(lambda z: C.phi_n(0, z), f, dx=1e-3)
        dEtaO = derivative(lambda z: C.eta_o(z), f, dx=1e-3)
    #     print 'dEtaO = ', dEtaO
        U = lambda z: C.U(z) + mn**4 * z**2 * (C.eta_s(z))/ (2* C.Cs)
        d2U = derivative(lambda z: U(z), f, dx=1e-3, n=2)
        d2EtaO = derivative(lambda z: C.eta_o(z), f, dx=1e-3, n=2)
        d2Phi = derivative(lambda z: C.phi_n(0, z), f, dx=1e-3, n=2)
    #     print n, d2Phi, d2EtaO, 2* C.Cs/mn**2 * d2Phi * C.phi_n(0,f)* I3num(n/2, meff)
        eta_contrib = derivative(lambda z: 1./C.eta_o(z), f, dx=1e-3, n=2)#2*dEtaO**2 / C.eta_o(f)**3 - d2EtaO / C.eta_o(f)**2
    #     print 'n= %f, f = %f, eta_contrib = %.12f'%(n, f, eta_contrib)
        res2 = -C.Cs / mn**4 #* C.phi_n(0, f) * dPhi/ (sqrt(pf**2 + meff**2) * mn**2)
        num = (C.Co*n*dEtaO/(mn**2 * C.eta_o(f)**2) - mn**2 * C.phi_n(0, f) * dPhi / sqrt(pf**2 + meff**2))**2

        res2 *= num
    #     res2 *= - mn**2 * C.phi_n(0, f) * dPhi / sqrt(pf**2 + meff**2)
        denom = (2 * C.Cs / mn**2 * dPhi**2 * I1num(n/2, meff) + d2U * C.Cs/mn**4 +
            C.Cs * C.Co * n**2 * eta_contrib / (2*mn**6) +
            2* C.Cs/mn**2 * d2Phi * C.phi_n(0,f) * I3num(n/2, meff))
        res2 /= denom
        mult = 1.
        if multiply:
            mult = 4 * pf * sqrt(pf**2 + meff**2) / (2 * pi**2)
        return f, (res + res2)*mult


    def f0(self, n, multiply=True):
        f = 0.
        res = []
        for _n in n:
            f, _res = self._f0(_n, f,multiply=multiply)
            res.append(_res)

        return np.array(res)

    def f0_der_NM(self, f, n, multiply = True):
        C = self.C
        mn = C.M[0]
        pf = eos.p_f(n)
        f, = eos.f_eq(np.array([n, 0.]), np.array([f]), 1, C)

        meff = mn*C.phi_n(0, f)

        print f
        mu = lambda z, f: derivative(lambda x: eos._E(np.array([f, x, 0.]), C), z, dx=1e-3)

        dmu_df = derivative(lambda x: mu(n, x), f, dx=1e-3)
        dmu_dn = derivative(lambda x: mu(x, f), n, dx=1e-3)
        d2E_df = derivative(lambda z: eos._E(np.array([z, n, 0.]), C), f, dx=1e-3, n=2)
    #     print 'dmu_dn = ', dmu_dn, ' ?= ', C.Co/ mn**2 /C.eta_o(f) + (3 * pi**2)**(2./3) * (n/2)**(-1./3) / sqrt(pf**2 + meff**2) / 6
    #     print derivative(lambda x: eos._E(np.array([f, x, 0.]), C), n, dx=1e-3, n=2)
        res = dmu_dn - (3 * pi**2)**(2./3) * (n)**(-1./3) / sqrt(pf**2 + meff**2) / 6 - dmu_df**2 / d2E_df

        mult = 4 * pf * sqrt(pf**2 + meff**2) / (2 * pi**2)
        if multiply:
            res *= mult
        return f, res

    def _f0_nm(self, n, f, multiply=True):
        C = self.C
        mn = C.M[0]
        pf = eos.p_f(n)
        
        f, = eos.f_eq(np.array([n, 0.]), np.array([f]), 1, C)
        print 'f=',f
        
        meff = mn * C.phi_n(0, f)
        
        Cv = lambda z: C.Co/C.eta_o(z) + C.Cr/(4*C.eta_r(z))
        res = Cv(f)/mn**2
        
        dPhi = derivative(lambda z: C.phi_n(0, z), f, dx=1e-3)
        
        dCv = derivative(lambda z: Cv(z), f, dx=1e-3)

        U = lambda z: C.U(z) + mn**4 * z**2 * (C.eta_s(z))/ (2* C.Cs)
        d2U = derivative(lambda z: U(z), f, dx=1e-3, n=2)
        d2Cv = derivative(lambda z: Cv(z), f, dx=1e-3, n=2)
        d2Phi = derivative(lambda z: C.phi_n(0, z), f, dx=1e-3, n=2)

#         den_func = lambda z: mn**4 * z**2 * C.eta_s(z) /2/C.Cs + C.U(z)
#         den_deriv = derivative(den_func, f, dx=1e-3, n=2) + d2Cv
        
        res2 = 1
        num = -(n*dCv/mn**2 - mn**2 * C.phi_n(0, f) * dPhi / sqrt(pf**2 + meff**2))**2
#         num = -meff**2/(pf**2 + meff**2) * C.Cs/mn**2
        res2 *= num
    
        denom = (mn**2 * dPhi**2 * I1num(n, meff) + d2U + d2Cv * n**2/2/mn**2 +
            mn**2 * d2Phi * C.phi_n(0,f) * I3num(n, meff))

#         denom = 1 + C.Cs/mn**2 * I1num(n, meff)
        res2 /= denom
        mult = 1.
        if multiply:
            mult = 2  * pf * sqrt(pf**2 + meff**2) / (2 * pi**2)
        return f, (res + res2)*mult
    
    def f0_nm(self, n, multiply=True):
        f = 0.
        res = []
        for _n in n:
            f, _res = self._f0_nm(_n, f, multiply=multiply)
            res.append(_res)

        return np.array(res)

    def _f1(self,n, f):
        f_f1 = f
        C = self.C
        pf = eos.p_f(n/2)
        mn = C.M[0]

        f_f1, = eos.f_eq(np.array([n/2, n/2]), np.array([f]), 1, C)
        f = f_f1
        meff = mn*C.phi_n(0,f)

        res = - C.Co * pf**2/(mn**2 * C.eta_o(f) * (pf**2 + meff**2))
        res /= (1. + 2 * C.Co/(mn**2 * C.eta_o(f)) * ( (2./3.) * I1num(n/2, meff) +
                                                       meff**2 * I4num(n/2, meff)))

        mult = 4 * pf * sqrt(pf**2 + meff**2) / (2 * pi**2)
        return f, res*mult

    def f1(self, n):
        f = 0.
        res = []
        for _n in n:
            f, f1 = self._f1(_n, f)
            res.append(f1)
        return np.array(res)

    def _f1_nm(self,n, f):
        f_f1 = f
        C = self.C
        pf = eos.p_f(n)
        mn = C.M[0]

        f_f1, = eos.f_eq(np.array([n, 0.]), np.array([f]), 1, C)
        f = f_f1
        meff = mn*C.phi_n(0,f)
        
        Cv = lambda z: C.Co/C.eta_o(z) + C.Cr/(4* C.eta_r(z))

        res = - Cv(f)/ mn**2 * pf**2/(pf**2 + meff**2)
        res /= (1. + Cv(f)/(mn**2) * ( (2./3.) * I1num(n, meff) +
                                                       meff**2 * I4num(n, meff)))

        mult = 2 * pf * sqrt(pf**2 + meff**2) / (2 * pi**2)
        return f, res*mult

    def f1_nm(self, n):
        f = 0.
        res = []
        for _n in n:
            f, f1 = self._f1_nm(_n, f)
            res.append(f1)
        return np.array(res)


    def UofE(self, i, n):
        f = eos.f_eq(n, np.array([self.C.f0]), 1, self.C)
        pots = eos.potentials(np.insert(n, 0, f), 5, self.C)
        print pots
        V = self.C.X_o[i] * pots[2]
        S =  (-self.C.M[i]) * (1 - self.C.phi_n(i, self.C.X_s[i]*self.C.M[0]/self.C.M[i]*f[0]))
        Uopt = lambda z: z + self.C.M[0] - sqrt((z + self.C.M[0] - V)**2 - S*(2 * self.C.M[i] + S))
        erange = np.linspace(-50./self.m_pi, 10., 100)
#         plt.plot(erange, self.m_pi*Uopt(erange))
#         plt.show()
        return self.m_pi*erange, self.m_pi*Uopt(erange)

    def dumpUofE(self, folderName, show=False):
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

        with open(join(folderName, 'UofE.dat'), 'w') as f:
            f.write(tab)
            
    def dumpUofE_anti(self, folderName, show=False):
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

        with open(join(folderName, 'UofE_anti.dat'), 'w') as f:
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

    def testHyper(self, z_l = 1., a_l = 1.):
        self.reset(hyper=1, nmin = 0., nmax = 8*self.n0, npoints=200)
        print self.E
        rho = self.concentrations()
        fig, ax = plt.subplots(1, 3)
        ax[0].plot(self.n/self.n0, rho)
        self.setDriver()
        n, m, r, mg = self.stars(npoints = 100)
        ax[1].plot(n/self.n0, m)
        ax2lines = []
        l, = ax[2].plot(self.n/self.n0, self.rho[:, 0])
        ax2lines.append(l)
        if self.C.sigma_kind == 1:
            l2, = ax[2].plot(self.n/self.n0, map(lambda z: self.C.Xs(2, z),
                                           self.rho[:,0]))
            ax2lines.append(l2)
            ax[2].legend(ax2lines, ['f(n)', 'x_sigma'], loc=0)

        plt.show()

    def diffArray(self, src, dx):
        res = []
        res.append((src[1] - src[0])/dx)
        _src = src[1:-1]
        for i, _s in enumerate(_src):
            res.append((src[i+2] - src[i])/(2*dx))
        n = src.size
        res.append((src[n-1]-src[n-2])/dx)
        return np.array(res)

    def P_chem(self, _n, leptons=True):
        sp = 1+self.C.sprime
        res = []
        for n in _n:
            mu_e = eos.mu(n, sp, self.C) - eos.mu(n, sp+1, self.C)
            n_e = 0.
            n_mu = 0.
            if (mu_e**2 - self.m_e**2 > 0):
                n_e += (mu_e**2 - self.m_e**2)**(1.5)/(3 * pi**2)

            if (mu_e**2 - self.m_mu**2 > 0):
                n_mu +=  (mu_e**2 - self.m_mu**2)**(1.5)/(3 * pi**2)

            if leptons:
                E = eos.E(n, self.C)
            else:
                E = eos._E(n, self.C)

            sum = 0.

            for i, r in enumerate(n[sp:]):
                sum += r *eos.mu(n, i+sp, self.C)

            if leptons:
                sum += n_e * mu_e + n_mu * mu_e

            res.append(sum - E)
        return np.array(res)

    def Jdiff(self, n):
        nf_n = np.array([n, 0.])
        nf_s = np.array([n/2, n/2])
        f_n = eos.f_eq(nf_n, np.array([self.C.f0]), 1, self.C)
        f_s = eos.f_eq(nf_s, np.array([self.C.f0]), 1, self.C)
#         print f_n, f_s
        EN = eos.EBind(np.concatenate((f_n, nf_n)), self.C)
        ES = eos.EBind(np.concatenate((f_s, nf_s)), self.C)
        return EN - ES

    def Ldiff(self):
        return 3*self.n0*derivative(lambda z: self.Jdiff(z),
                                    self.n0, dx=1e-3)

    def Ksymm_diff(self):
        return 9*self.n0**2 * derivative(lambda z: self.Jdiff(z),
                                         self.n0, dx=1e-3, n=2)

    def dumpChi(self, folderName):
        scalings = self.dumpScalings(None)
        chi_s = []
        chi_o = []
        chi_r = []

        for i, _n in enumerate(scalings[:,0]):
            chi_s.append(self.C.phi_n(0, scalings[i, 1]) / sqrt(scalings[i, 2]))
            chi_o.append(self.C.phi_n(0, scalings[i, 1]) / sqrt(scalings[i, 3]))
            chi_r.append(self.C.phi_n(0, scalings[i, 1]) / sqrt(scalings[i, 4]))

        tab = np.array([scalings[:, 0], chi_s, chi_o, chi_r]).transpose()
        table = tabulate(tab, ['n/n0', 'chi_s', 'chi_o', 'chi_r'],
                         tablefmt='plain')

        with open(join(folderName, 'chi_N.dat'), 'w') as f:
            f.write(table)

        return tab

    def dumpStarProfile(self, n_c, folderName=None, file_suffix=None):
        n, m, r, mg1, mg2 = self.stars_crust(nmin=2*self.n0, nmax=n_c,
                                             npoints=2)

        dr = self.dr
        nSize = self.dr.nSize
        N = dr.getLastN(nSize)[:-1]
        R = dr.getLastR(nSize)[:-1]
        M = dr.getLastM(nSize)[:-1]
        E = dr.getLastE(nSize)[:-1]
        P = dr.getLastP(nSize)[:-1]

        if folderName is not None:
            if not os.path.exists(folderName):
                os.makedirs(folderName)

            fname = 'starDump_n=%.2f'%(n_c/self.n0)
            if file_suffix is not None:
                fname+=file_suffix

            fname += '.dat'
            tab = np.array([N/self.n0, M, R, E, P]).transpose()
            table = tabulate(tab, ['n/n_0', 'M/M_sun', 'R [km]', 'E [m_pi^4]',
                                   'P [m_pi^4]'], tablefmt='plain')

            with open(join(folderName, fname), 'w') as f:
                f.write(table)

        return N, M, R, E, P
    
    def dumpJ(self, folderName):
        n = np.linspace(0., 8*self.n0, 800, endpoint=0)
        f = 0.
        res = []
        for _n in n:
            f, = eos.f_eq(np.array([_n/2, _n/2]), np.array([f]), 1, self.C)
            res.append(eos.J(_n, self.C, f))
        res = np.array(res)
        tab = np.array([n/self.n0, res]).transpose()
        table = tabulate(tab, ['n/n_0', 'J'], tablefmt='plain')
        with open(join(folderName, 'J.dat'), 'w') as f:
            f.write(table)
        return n, res

    def dumpVs(self, folderName=None):
#         n = np.linspace(0., 8*self.n0, 80, endpoint=0)
#         Es = self.Esymm(n)
#         Ps = self.Psymm(n)/self.const/self.m_pi**4
#         Vs = np.diff(Ps)/np.diff(Es)
        
        self.reset(hyper=0, nmin=0., nmax=8*self.n0)
        
        vsNs = np.gradient(self.P)/np.gradient(self.E)
        tab = np.array([self.n[:]/self.n0, vsNs]).transpose()
        table = tabulate(tab, ['n/n0', 'vsNs'], tablefmt='plain')
        
        if folderName is not None:
            with open(join(folderName, 'vsNs.dat'), 'w') as f:
                f.write(table)
#             print self.n, type(self.n)
            es = self.Esymm(self.n)
            ps = self.Psymm(self.n)/self.const/self.m_pi**4
            vs = np.gradient(ps)/np.gradient(es)
            tab = np.array([self.n[:]/self.n0, vs]).transpose()
            table = tabulate(tab, ['n/n0', 'vs_sym'], tablefmt='plain')
            with open(join(folderName, 'vs_sym.dat'), 'w') as f:
                f.write(table)
            
                    
        else:
            plt.plot(self.n[1:]/self.n0, vsNs)
            plt.ylim([0., 1])
            plt.show()  
            
    def dumpVsHyper(self, folderName=None, show=0, npoints=1000):
#         n = np.linspace(0., 8*self.n0, 80, endpoint=0)
#         Es = self.Esymm(n)
#         Ps = self.Psymm(n)/self.const/self.m_pi**4
#         Vs = np.diff(Ps)/np.diff(Es)
        self.C.phi_meson = 0
        self.reset(hyper=1, nmin=0., nmax=8*self.n0, npoints=npoints, timeout=6)
        
        vsNs = np.gradient(self.P)/np.gradient(self.E)
        tab = np.array([self.n[:]/self.n0, vsNs]).transpose()
        table = tabulate(tab, ['n/n0', 'vs^2'], tablefmt='plain')
        
        if folderName is not None:
            with open(join(folderName, 'vsHyper.dat'), 'w') as f:
                f.write(table)            
                    
        if show:
            plt.plot(self.n[1:]/self.n0, vsNs)
            plt.ylim([0., 1])
            plt.show()
            
        self.C.phi_meson = 1
        self.reset(hyper=1, nmin=0., nmax=8*self.n0, npoints=npoints, timeout=6)
        
        vsNs = np.gradient(self.P)/np.gradient(self.E)
        tab = np.array([self.n[:]/self.n0, vsNs]).transpose()
        table = tabulate(tab, ['n/n0', 'vs^2'], tablefmt='plain')
        
        if folderName is not None:
            with open(join(folderName, 'vsHyperPhi.dat'), 'w') as f:
                f.write(table)             
        
        if show:
            plt.plot(self.n[1:]/self.n0, vsNs)
            plt.ylim([0., 1])
            plt.show()
        
        self.C.phi_meson = 0
        
            
    def dumpVsSymLow(self, folderName=None):
        n = np.linspace(0., 0.01, 1000)
        es = self.Esymm(n)
        ps = self.Psymm(n)/self.const/self.m_pi**4
        vs = np.gradient(ps)/np.gradient(es)
#         plt.plot(n[1:]/self.n0, vs)
#         plt.show() 
        tab = np.array([n[:]/self.n0, vs]).transpose()
        table = tabulate(tab, ['n/n0','vs'], tablefmt='plain')
        if folderName is not None:
            with open(join(folderName, 'vs_sym_detail.dat'), 'w') as f:
                f.write(table)
        else:
            plt.plot(n[1:]/self.n0, vs)
            plt.show()
        return n[1:], vs
            
    def getContrib(self):
        om = []
        rho = []
#         xo = np.array([self.C.X_o[i] for i in range(8)])
#         xr = np.array([self.C.X_o[i] for i in range(8)])
#         t3 = np.array([self.C.T[i] for i in range(8)])
        for r in self.rho:
            om_sum = 0.
            rho_sum = 0.
            f = r[0]
            for i, _r in enumerate(r[1:]):
                om_sum += self.C.X_o[i] * _r
                rho_sum += self.C.X_r[i] * _r * self.C.T[i]
            om.append(self.C.Co*om_sum**2 / (2 * self.C.M[0]**2 * self.C.eta_o(f)))
            rho.append(self.C.Cr*rho_sum**2 / (2 * self.C.M[0]**2 * self.C.eta_r(f)))
        return np.array([om, rho]).transpose()
    
    def eval_DU(self, crust=0):
        if crust:
            n,m,r,mg,mg2 = self.stars_crust()
        else:
            self.reset(npoints = 400)
            self.setDriver()
            n,m,r,mg = self.stars(npoints=50)
        np = self.concentrations()[:, 1]
        im = interp1d(n, m)
        inp = interp1d(self.n, np)
        
        x_du=[]
        l_n_e = []
        l_n_mu = []
        for _mu in self.mu_e:
            n_e = 0
            n_mu = 0
            if _mu > self.m_e:
                n_e += (_mu**2 - self.m_e**2)**(1.5)/(3*pi**2)
            if _mu > self.m_mu:
                n_mu += (_mu**2 - self.m_mu**2)**(1.5)/(3*pi**2)
            
            l_n_e.append(n_e)
            l_n_mu.append(n_mu)
            
            if n_e + n_mu > 0:
                x_du.append(1./(1 + (1 + (n_e/(n_e + n_mu))**(1./3))**3))
            else:
                x_du.append(0.11)
        ix_du = interp1d(self.n, x_du)
        res = optimize.minimize(lambda z: (inp(z) - ix_du(z))**2, [2.5*self.n0])
        x = res.x[0]

        print res
        print x, x/self.n0
        print im(x)
        plt.plot(self.n, inp(self.n), self.n, ix_du(self.n))
        plt.show()
        return [x/self.n0, im(x)]
    
    def dumpHyperScalings(self, folderName):
        if not os.path.exists(folderName):
            os.makedirs(folderName)
            
        frange = np.linspace(0, 1, 100)
        phi_sc = map(self.C.eta_p, frange)
        chi_p = sqrt((1 - frange)**2/np.array(phi_sc))
        etap_tab = np.array([frange, phi_sc]).transpose()
        etap_table = tabulate(etap_tab, ['f', 'eta_p'], tablefmt='plain')
        with open(join(folderName, 'etap.dat'), 'w') as f:
            f.write(etap_table)
#         print chi_p
        chi_tab = np.array([frange, chi_p]).transpose()
#         print chi_tab
        chi_table = tabulate(chi_tab, ['f', 'chi_p'], tablefmt='plain')
        with open(join(folderName, 'chi_p.dat'), 'w') as f:
            f.write(chi_table)
        
        tab = []
        for _f in frange:
            _tab = []
            _tab.append(_f)
            for i in xrange(2, 8):
                _tab.append(self.C.Xs(i, _f)/self.C.X_s[i])
            tab.append(_tab)
        etap_tab = np.array(tab)#.transpose()
        etap_table = tabulate(etap_tab, ['f', 'L','S-','S0','S+','X-','X0'], tablefmt='plain')
        with open(join(folderName, 'xi_s.dat'), 'w') as f:
            f.write(etap_table)
            
        self.reset(hyper=1, npoints=200, timeout=6)
        frange = self.rho[:, 0]
        phi_sc = map(self.C.eta_p, frange)
        etap_tab = np.array([self.n/self.n0, phi_sc]).transpose()
        etap_table = tabulate(etap_tab, ['n/n0', 'eta_p'], tablefmt='plain')
        with open(join(folderName, 'etap_N.dat'), 'w') as f:
            f.write(etap_table)
        
        chi_p = sqrt((1 - frange)**2/np.array(phi_sc))
        chi_tab = np.array([self.n/self.n0, chi_p]).transpose()
        chi_table = tabulate(chi_tab, ['n/n0', 'chi_p'], tablefmt='plain')
        with open(join(folderName, 'chi_p_N.dat'), 'w') as f:
            f.write(chi_table)
        
        tab = []
        for j, _f in enumerate(frange):
            _tab = []
            _tab.append(self.n[j]/self.n0)
            for i in xrange(2, 8):
                _tab.append(self.C.Xs(i, _f)/self.C.X_s[i])
            tab.append(_tab)
        etap_tab = np.array(tab)#.transpose()
        etap_table = tabulate(etap_tab, ['n/n0', 'L','S-','S0','S+','X-','X0'], tablefmt='plain')
        with open(join(folderName, 'xi_s_N.dat'), 'w') as f:
            f.write(etap_table)
            
        frange = np.linspace(0, 1, 100)
        self.C.phi_meson=1
        
        self.reset(hyper=1, npoints=400, timeout=6)
        frange = self.rho[:, 0]
        phi_sc = map(self.C.eta_p, frange)
        etap_tab = np.array([self.n/self.n0, phi_sc]).transpose()
        etap_table = tabulate(etap_tab, ['n/n0', 'eta_p'], tablefmt='plain')
        with open(join(folderName, 'etap_N_phi.dat'), 'w') as f:
            f.write(etap_table)
        
        chi_p = sqrt((1 - frange)**2/np.array(phi_sc))
        chi_tab = np.array([self.n/self.n0, chi_p]).transpose()
        chi_table = tabulate(chi_tab, ['n/n0', 'chi_p'], tablefmt='plain')
        with open(join(folderName, 'chi_p_N.dat'), 'w') as f:
            f.write(chi_table)
        
        tab = []
        for j, _f in enumerate(frange):
            _tab = []
            _tab.append(self.n[j]/self.n0)
            for i in xrange(2, 8):
                _tab.append(self.C.Xs(i, _f)/self.C.X_s[i])
            tab.append(_tab)
        etap_tab = np.array(tab)#.transpose()
        etap_table = tabulate(etap_tab, ['n/n0', 'L','S-','S0','S+','X-','X0'], tablefmt='plain')
        with open(join(folderName, 'xi_s_N_phi.dat'), 'w') as f:
            f.write(etap_table)   
        
    
class EosConstructor(object):
    def __init__(self):
        self.Fs = None
        self.Fn = None
        self.Fx = None
        self.eta_o = lambda z: 1
        self.eta_r = lambda z: 1
        self.setConstants()
        self.setParams()

    def setFs(self, Fs, n=None):
        if isinstance(Fs, interp1d):
            self.Fs = Fs
        else:
            if n is not None:
                self.Fs = interp1d(n, Fs, kind='cubic')
            else:
                raise RuntimeError("""setFs requires either spline or x
                and y data""")

    def getFs(self):
        assert self.Fs is not None
        return self.Fs

    def setFn(self, Fn, n=None):
        if isinstance(Fn, interp1d):
            self.Fn = Fn
        else:
            if n is not None:
                self.Fn = interp1d(n, Fn, kind='cubic')
            else:
                raise RuntimeError("""setFn requires either spline or x
                        and y data""")

    def getFn(self):
        assert self.Fn is not None
        return self.Fn

    def setScaling(self, eta_o, eta_r):
        assert isinstance(eta_o, LambdaType) and isinstance(eta_r, LambdaType)
        self.eta_o = eta_o
        self.eta_r = self.eta_r

    def getFx(self):
        pass

    def setConstants(self, Co=54.6041, Cr=121.69):
        self.Co = Co
        self.Cr = Cr

    def setParams(self, mn=938./135, m_pi=135.):
        self.mn=mn
        self.m_pi = m_pi

    def checkReady(self):
        assert self.Fn is not None
        assert self.Fs is not None
        assert isinstance(self.eta_o, LambdaType)
        assert isinstance(self.eta_r, LambdaType)

    def pf(self, n):
        return (3*pi**2*n)**(1./3.)

    def ES(self, n):
        f = self.Fs
        C = lambda z: self.Co / self.eta_o(f(z))
        res = integrate.quad(lambda z: sqrt(self.pf(z/2)**2 +
                                            self.mn**2*(1-f(z))**2) +
                             C(z) * z/self.mn**2, 0., n, epsrel=1e-11)[0]
        return res

    def EN(self, n):
        f = self.Fn
        C = lambda z: self.Co / self.eta_o(f(z))+ self.Cr / 4. / self.eta_r(f(z))
        res = integrate.quad(lambda z: sqrt(self.pf(z)**2 +
                                            self.mn**2*(1-f(z))**2) +
                             C(z) * z/self.mn**2, 0., n, epsrel=1e-11)[0]
        return res

    def PN(self, n):
        return n*derivative(lambda z: self.EN(z), n, dx=1e-4) - self.EN(n)

class FromEos(EosConstructor):
    def __init__(self):
        super(FromEos, self).__init__()

    def FsFromEos(self, E, P, N):
        assert E.shape == P.shape == N.shape
        res = []
        fEs = 0.
        C = lambda z: self.Co / self.eta_o(z)
        for i, _n in enumerate(N):
            fEs = root(lambda z: z - (1 - sqrt((E[i] + P[i] - C(z)
             * _n**2 / self.mn**2)**2 / (_n**2) - self.pf(_n/2)**2)/self.mn), [fEs], tol=1e-35).x
            res.append(fEs)
        res = np.array(res)[:,0]
        print res
        assert res.shape == N.shape
        self.Fs = interp1d(N, res, kind='cubic')
        return res

    def FnFromEos(self, E, P, N):
        assert E.shape == P.shape == N.shape
        res = []
        fEn = 0.
        C = lambda z: self.Co / self.eta_o(z) + self.Cr / 4 / self.eta_r(z)
        for i, _n in enumerate(N):
            fEn = root(lambda z: z - (1 - sqrt((E[i] + P[i] - C(z)
             * _n**2 / self.mn**2)**2 / (_n**2) - self.pf(_n)**2)/self.mn), [fEn], tol=1e-21).x
            res.append(fEn)

        res = np.array(res)[:,0]
        self.Fn = interp1d(N, res, kind='cubic')
        return res
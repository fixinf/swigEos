from scipy.interpolate.interpolate import interp1d
from scipy.misc.common import derivative
from scipy.optimize.optimize import fmin, fmin_cg, fmin_powell, fmin_bfgs
from scipy.optimize.zeros import brentq
from veusz.utils.utilfuncs import isiternostr
import eosWrap as eos

__author__ = 'const'
import numpy as N
from scipy.optimize import root
from math import pi, asinh, sqrt, tanh
from matplotlib import pyplot as plt

import lmfit as lm

class NLVector():
    def __init__(self, Cs, Co, Cr, b, c, Co4, Cr4, cut_a=0., cut_b=70., cut_c=100.):
        self.Cs = Cs
        self.Co = Co
        self.Cr = Cr
        self.b = b
        self.c = c
        self.Co4 = Co4
        self.Cr4 = Cr4
        self.cut_a = cut_a
        self.cut_b = cut_b
        self.cut_c = cut_c
        self.mpi = 135
        self.mpi4_2_mevfm3 = 135 * (135 / 197.33)**3
        self.mn = 939/self.mpi
        self.mo = 782.501/self.mpi
        self.mr = 763/self.mpi
        self.go = sqrt(self.Co) / self.mn * self.mo
        self.gr = sqrt(self.Cr) / self.mn * self.mr
        # self.n0 = ((1.3*197.33/self.mpi)**3) / ((3*np.pi**2)/2.)
        self.n0 = 0.148 * (197.33 / 135)**3
        self.f_last = 0.2
        self.f0 = self.f_eq(self.n0/2, self.n0/2)
        self.m_e = 0.5/self.mpi
        self.m_mu = 105.0 / self.mpi
        self.vec_last = [0, 0]
        self.np_last = 0.

    # def __init__(self):
    #     self.mpi = 135
    #     self.mn = 939/self.mpi
    #     self.mo = 782.5/self.mpi
    #     self.mr = 763/self.mpi
    #     self.n0 = ((1.42*197.33/self.mpi)**3) / ((3*np.pi**2)/2.)

    def pf(self, n):
        if n < 0:
            return 0
        return (3 * pi**2 * n)**(1/3)

    def U(self, f):
        fs = 100500.
        if hasattr(self, 'f0'):
            fs = self.f0 + self.cut_c * (1 - self.f0)
        return self.mn**4 * (self.b * f**3 /3 + self.c * f**4 / 4) + \
                self.cut_a/2 * (1 + tanh(self.cut_b * (f - fs)))


    def I1(self, m, n):
        pf = self.pf(n)
        return (-m**4*asinh(pf/m)/8 + m**3*pf/(8*sqrt(1 + pf**2/m**2))
                + 3*m*pf**3/(8*sqrt(1 + pf**2/m**2)) +
                pf**5/(4*m*sqrt(1 + pf**2/m**2)))/pi**2


    def I2(self, m, n):
        pf = self.pf(n)
        return (-m**3*asinh(pf/m)/2 + m**2*pf/(2*sqrt(1 + pf**2/m**2)) +
                   3*pf**3/(8*sqrt(1 + pf**2/m**2)) + pf**3/(8*(1 + pf**2/m**2)**(3/2)) -
                   pf**5/(4*m**2*sqrt(1 + pf**2/m**2)) + 3*pf**5/(8*m**2*(1 + pf**2/m**2)**(3/2)) +
                   pf**7/(4*m**4*(1 + pf**2/m**2)**(3/2)))/pi**2



    def f_eq(self, nn, np):
        meff = lambda x: self.mn * (1 - x)
        eq = lambda z: self.mn**4 * z / self.Cs - self.mn * self.I2(meff(z), nn) - self.mn * self.I2(meff(z), np)\
                       + derivative(lambda x: self.U(x), z, dx=1e-3)
        # frange = N.linspace(0, 1, 100)
        # plt.plot(frange, list(map(eq, frange)))
        # plt.show()
        self.f_last = root(eq, self.f_last).x[0]
        return self.f_last

    def vec_eq(self, nn, np):
        if isiternostr(nn):
            nn = nn[0]
        if isiternostr(np):
            np = np[0]
        eq = lambda x: [-self.mo**2 * x[0] + self.go * (nn + np) - 4*self.Co4 * x[0]**3 - 2 * self.Cr4 * x[1]**2 * x[0],
                            -self.mr**2 * x[1] + self.gr * (np - nn)/2 - 2 * self.Cr4 * x[1] * x[0]**2]
        self.vec_last = root(eq, self.vec_last).x
        return self.vec_last

    def _E(self, f, om, rho, nn, np):
        E = self.mn**4 * f**2 / (2 * self.Cs)
        E += -self.mo**2 * om**2 /2 + self.go * om * (nn + np) - self.Co4 * om**4
        E += -self.mr**2 * rho**2 /2 + self.gr * rho * (np - nn)/2 - self.Cr4 * om**2 * rho**2
        E += self.U(f)
        E += self.I1(self.mn*(1-f), nn)
        E += self.I1(self.mn*(1-f), np)
        return E

    def E(self, nn, np):
        f = self.f_eq(nn, np)
        om, rho = self.vec_eq(nn, np)
        return self._E(f, om, rho, nn, np)

    def mu(self, nn, np):
        # if nn < 0 or np < 0:
        #     return [0.0, 0.0]
        return [derivative(lambda z: self.E(z, np), nn, dx=1e-5, order=7),
                derivative(lambda z: self.E(nn, z), np, dx=1e-5, order=7)]

    def np_eq(self, n):
        def fun(params, n):
            np = params['np'].value
            mu = self.mu(n - np, np)
            mu_e = mu[0] - mu[1]
            n_e = 0.
            n_mu = 0.
            n_l = 0.
            if mu_e**2 - self.m_e**2 > 0:
                n_e += (mu_e**2 - self.m_e**2)**(1.5)/(3 * pi**2)
            if mu_e**2 - self.m_mu**2 > 0:
                n_mu += (mu_e**2 - self.m_mu**2)**(1.5)/(3 * pi**2)
            # if mu_e > 0:
            n_l = n_e + n_mu
            # print('fun:', n, np, n_l)
            return [(np - n_l), n_e, n_mu]

        fit_params = lm.Parameters()
        fit_params.add('np', value=self.np_last, min=0, max=n)
        min = lm.minimize(lambda z: N.array([fun(z, n)[0]]), fit_params, method='leastsq')#, **{'xtol': 1e-16, 'ftol': 1e-16})
        # print(min.message)
        # lm.report_fit(fit_params)
        self.np_last = fit_params['np'].value
        # try:
        #     # res = brentq(lambda z: fun(n - z, z)[0], 0., n)#, 0.1, tol=1e-8, options={'epsfcn' : 1e-4}).x[0]
        #     self.np_last = fmin(lambda z: fun(n - z, z)[0], self.np_last, xtol=1e-8)#, options={'epsfcn' : 1e-4}).x[0]
        # except ValueError:
        #     res = 0.
        return [self.np_last] + fun(fit_params, n)[1:]

    def K(self, n):
        return self.mpi * 9*n**2 * derivative(lambda z: self.E(z/2, z/2)/z, n,
                                              dx=1e-3, n=2, order=9)

    def J(self, n):
        return self.mpi * n / 8 * derivative(lambda z: self.E(n/2 - z, n/2 + z), 0.,
                                         dx=1e-3, n=2, order=9)

    def P(self, nrange, E):
        return nrange * N.gradient(E, [nrange[1]-nrange[0]]) - E

    def EPN_NS(self, npoints, ret_np=0):
        res = []
        nplist = []
        nrange = N.linspace(0, 10*self.n0, npoints, endpoint=1)
        for n in nrange:
            np, ne, nmu = self.np_eq(n)
            print(np, ne, nmu)
            nplist.append(np/n)
            e = self.E(n-np, np) + self.I1(self.m_e, ne) + self.I1(self.m_mu, nmu)
            res.append(e)
        E = N.array(res)
        P = self.P(nrange, E)
        if not ret_np:
            return E, P, nrange
        else:
            return E, P, nrange, N.nan_to_num(N.array(nplist))

    def tabEos(self, name='test', npoints=80):
        E, P, n, np = self.EPN_NS(npoints, ret_np=1)
        nsym = N.array([[z/2, z/2] for z in n])
        Esym = N.nan_to_num(N.array([self.E(*z) for z in nsym]))
        Ebind = self.mpi * (Esym/n - self.mn)
        E *= self.mpi4_2_mevfm3
        P *= self.mpi4_2_mevfm3
        n /= self.n0

        N.savetxt(name+'_eos.dat', N.array([n, E, P, np, Ebind]).transpose(), fmt='%.6e')
        return [E, P, n]


    def tabMass(self, name='test', nstars=40):
        dr = eos.KVDriver()
        E, P, n = self.tabEos(name=name)
        n_dense = N.linspace(0, 10*self.n0, 800)
        kind = 'cubic'
        E_dense = interp1d(n*self.n0, E, kind=kind)(n_dense)
        P_dense = interp1d(n*self.n0, P, kind=kind)(n_dense)
        dr.set(E_dense/self.mpi4_2_mevfm3, P_dense/self.mpi4_2_mevfm3, n_dense)
        stars = N.linspace(0, 10*self.n0, nstars)
        MR = []
        for n in stars:
            MR.append(eos.star_crust2(n, 3, dr, 1e-11))
        MR = N.array(MR)
        N.savetxt(name+'_mass.dat', N.insert(MR, 0, stars/self.n0, axis=1), fmt='%.4f')
        return MR


    def tabMassCrust(self, name='test', nstars=100, ncut_crust=.45, ncut_eos=.6):
        dr = eos.KVDriver()
        E, P, n = self.tabEos(name=name)
        n_dense = N.linspace(0, 10*self.n0, 10000)
        kind = 'linear'

        e = []
        p = []
        ncrust = []
        with open("/home/const/workspace/swigEosWrapper/crust.dat", 'r') as f:
            for line in f:
                # print line
                _e, _p, _n = line.split()
                if float(_n) < ncut_crust:
                    e.append(float(_e)/self.mpi**4 * self.mpi4_2_mevfm3)
                    p.append(float(_p)/self.mpi**4 * self.mpi4_2_mevfm3)
                    ncrust.append(float(_n)*self.n0)
        e = N.array(e)
        p = N.array(p)
        ncrust = N.array(ncrust)
        print(ncrust)

        # gamma = 1. / 4.
        # n_dense = N.linspace(0., (10*self.n0) ** gamma, 1000)
        # n_dense = n_dense ** (1. / gamma)

        k_eos = N.where(n > ncut_eos)[0][0]

        Enew = N.concatenate((e, E[k_eos:]))
        Pnew = N.concatenate((p, P[k_eos:]))
        Nnew = N.concatenate((ncrust, n[k_eos:]*self.n0))

        # plt.plot(n, P, label='eos')
        # plt.plot(ncrust/self.n0, p, label='crust')
        # plt.plot(Nnew/self.n0, Pnew, label='tot')
        # plt.legend()
        # plt.show()

        E_dense = interp1d(Nnew, Enew, kind=kind)(n_dense)
        P_dense = interp1d(Nnew, Pnew, kind=kind)(n_dense)

        # plt.plot(n_dense/self.n0, P_dense)
        # plt.show()
        dr.set(E_dense/self.mpi4_2_mevfm3, P_dense/self.mpi4_2_mevfm3, n_dense)
        stars = N.linspace(.4*self.n0, 10*self.n0, nstars)
        MR = []
        for n in stars:
            MR.append(eos.star_crust2(n, 3, dr, 1e-11))
        MR = N.array(MR)
        N.savetxt(name+'_mass_crust.dat', N.insert(MR, 0, stars/self.n0, axis=1), fmt='%.4f')
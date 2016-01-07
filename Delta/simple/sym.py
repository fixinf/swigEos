# Nonlinear Walecka model Scipy realisation with Delta (1232) isobars
#
from scipy.constants.constants import pi
from math import pi, asinh, sqrt, tanh
from scipy.misc.common import derivative
from scipy.optimize import optimize
from scipy.optimize._root import root
from scipy.optimize.minpack import leastsq

# import eosWrap as eos

__author__ = 'const'
import numpy as np
from matplotlib import pyplot as plt
import matplotlib


matplotlib.rcParams.update({'font.size': 22})
matplotlib.rcParams.update({'lines.linewidth': 3})
mpi = 135.
n0 = 0.17 * (197.33/mpi)**3
mn = 939./mpi
md = 1232/mpi

Cs = 246
Co = 156.3
b = 1.8e-3
c = 2.87e-4

xo = 1
xs = 1.3

def p_f(n, s):
    if n > 0:
        return (6 * pi**2 * n/ (2*s+1))**(1/3)
    else:
        return 0

def I1( m, n, s):
    pf = p_f(n, s)
    return (-m**4*asinh(pf/m)/8 + m**3*pf/(8*sqrt(1 + pf**2/m**2))
            + 3*m*pf**3/(8*sqrt(1 + pf**2/m**2)) +
            pf**5/(4*m*sqrt(1 + pf**2/m**2)))/pi**2

def U(f):
    return mn**4*(b*f**3 / 3 + c*f**4/4)

def E(nn, nd, f):
    Es = mn**4 * f**2 / (2 * Cs)
    Eo = Co * (nn + xo*nd)**2 / (2*mn**2)
    Eu = U(f)
    Ekn = 2*I1(mn*(1-f), nn/2, .5)
    Ekd = 8*I1(md*(1-mn/md * xs * f), nd/4, 1.5)
    E = Es + Eo + Eu + Ekn + Ekd
    # print(Es, Eo, Eu, Ekn, Ekd)
    return E

def f_eq(nn, nd, finit=.1):
    def eq(nn, nd, f):
        return derivative(lambda z: E(nn, nd, z), f, dx=1e-4, order=7)
    # return root(lambda z: eq(nn, nd, z), x0=finit, tol=1e-12)['x']
    return leastsq(lambda z: eq(nn, nd, z), x0=finit, xtol=1e-12)[0]

def Eb_N(n):
    e = []
    f = []
    _f = 0.
    for _n in n:
        _f = f_eq(_n, 0, finit=_f)[0]
        f.append(_f)
        e.append(E(_n, 0, _f))
    return mpi*(np.array(e)/n - mn), np.array(f)

def mu_d(n, nd):
    nn = n - nd
    f = f_eq(nn, nd, finit=.5)[0]
    return derivative(lambda z: E(n - nd, z, f), nd)

def mu_n(n, nd):
    nn = n - nd
    f = f_eq(nn, nd, finit=.3)[0]
    return derivative(lambda z: E(z, nd, f), nn, dx=1e-3)

def nd_eq(n, ninit=0.0, finit=0.1):
    def eq(nn, nd):
        if nn < 0:
            return [-100500]
            # pass
        f = f_eq(nn, nd, finit=finit)[0]
        res = (derivative(lambda z: E(z, nd, f), nn, dx=1e-3, order=3) -
                derivative(lambda z: E(nn, z, f), nd, dx=1e-3, order=3))
        # print(nn, nd, f, res)
        # print([res, f])
        return [res, f, derivative(lambda z: E(nn, z, f), nd, dx=1e-5, order=7)]
    zrange = np.linspace(0, n, 100)

    # print(eq(n/2, n/2))
    res = leastsq(lambda z: eq(n-z, z)[0], x0=ninit)[0][0]
    # if res > 1e-2:
        # plt.plot(zrange, list(map(lambda z: eq(n-z, z)[0], zrange)))
        # plt.show()
    # if res < 1e-6:
    #     res = 0.
    return res, eq(n-res, res)[1], eq(n-res, res)[0], eq(n-res, res)[2]

def get_nd(nrange):
    nd = []
    mu = []
    _nd = 0.
    _f = 0.
    for n in nrange:
        res = nd_eq(n, ninit=_nd, finit=_f)
        print(res)
        _nd = res[0]
        _f = res[1]
        mu.append(res[3])
        if _nd < 0:
            _nd = 0.
        nd.append(_nd)
    return np.array(nd), np.array(mu)

def Eb(nrange, nd):
    f = []
    _e = []
    _f = 0
    for i,_n  in enumerate(nrange):
        _f = f_eq(_n-nd[i], nd[i], finit=_f)[0]
        f.append(_f)
        _e.append(E(_n-nd[i], nd[i], _f))
    return mpi*(np.array(_e)/nrange - mn), np.array(f)

print(mu_n(n0, 0.)*135)
f0 = f_eq(n0, 0)
print(f0)
pf0 = p_f(n0, .5)
Ef0 = sqrt(pf0**2 + mn**2 * (1-f0)**2) + Co * n0 / mn**2
print(Ef0*mpi)
exit()


print(f_eq(2.69696969697, 0.3, finit=.6))
# exit()
nrange = np.linspace(0*n0, 18*n0, 200, endpoint=0)
ndd = np.array([0. for n in nrange])
# e,f = Eb(nrange, ndd)
# plt.plot(nrange/n0, e)
# plt.xlim([0.6, 3.6])
# plt.ylim([-30, 80])
# plt.show()
# exit()
# print(nd_eq(3., finit=.6))
ndd, mu = get_nd(nrange)
plt.plot(nrange/n0, [0.857 for n in nrange])
plt.plot(nrange/n0, ndd/nrange)
plt.show()
# exit()
e,f = Eb(nrange, ndd)
# plt.plot(nrange/n0, e)
# plt.xlim([0.6, 3.6])
# plt.ylim([-30, 80])
# plt.show()
#
_n = n0

print(E(_n, _n, .6))
#### compare fermi momenta
print('p_f:')
s = 1.5
print(p_f(n0, s))
print(eos.p_f(n0, 2*s+1))

#### compare kinetic integrals
print('E_kin:')
print(8*I1(mn, n0/4, s))
print(4*eos.kineticInt(n0/4, mn, 2*s+1))

####


# exit()
np.savetxt("out_delta.dat", np.array([nrange/n0, e, f, ndd/nrange, mpi*mu]).transpose(),fmt='%.6f')
mu_dd = [derivative(lambda z: mpi*E(n, z, f_eq(n, n)), n, dx=1e-3) for n in nrange/2]
np.savetxt("mu_test.dat", np.array([nrange/n0, mu_dd]).transpose(), fmt='%.6f')
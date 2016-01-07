from math import sqrt
from numpy.core.umath import sign
from scipy import optimize
from scipy.misc.common import derivative

__author__ = 'const'
import Models2
import eosWrap as eos
import numpy as np
from matplotlib import pyplot as plt

m = Models2.KVOR()
C = m.C


def mu_func(i, n_in, rho_0, rho_c, mu_c, z):
    n = np.copy(n_in)
    n[i] = z
    return eos.E_rho(n, rho_0, rho_c, mu_c[0], C)


def mu_deriv(i, n_in, rho_0, rho_c, mu_c):
    return derivative(lambda z: mu_func(i, n_in, rho_0, rho_c, mu_c, z), n_in[i], dx=1e-4, order=5)


def mu_ch(n_in, rho_0, rho_c):
    def eq(n_in, rho_0, rho_c, z):
        mu_e = mu_deriv(1, n_in, rho_0, rho_c, z) - mu_deriv(2, n_in, rho_0, rho_c, z)
        return mu_e - z

    return optimize.leastsq(lambda z: eq(n_in, rho_0, rho_c, z), [0.])[0]

def rho_eq(n_n, n_p):
    def eq(n_n, n_p, rho):
        rho_0 = rho[0]
        rho_c = rho[1]
        f = f_eq(n_n, n_p, rho_0, rho_c)[0]
        n_in = np.array([f, n_n, n_p])
        mu_r = mu_ch(n_in, rho_0, rho_c)[0]
        return [derivative(lambda z: eos.E_rho(n_in, z, rho_c, mu_r, C), rho_0, dx=1e-4, order=3),
                derivative(lambda z: eos.E_rho(n_in, rho_0, z, mu_r, C), rho_c, dx=1e-4, order=3)]
        # return [derivative(lambda z: eos.E_rho(n_in, z, 0., C), rho_0, dx=1e-4, order=7)]
    return optimize.root(lambda z: eq(n_n, n_p, z), [0.2, 1.])['x']

def f_eq(n_n, n_p, rho_0, rho_c):
    n_in = np.array([n_n, n_p, 0.])
    def eq(f):
        f = f[0]
        n_in = np.array([f, n_n, n_p])
        mu_c = mu_ch(n_in, rho_0, rho_c)[0]
        return derivative(lambda z: eos.E_rho(np.array([z, n_n, n_p]), rho_0, rho_c, mu_c, C), f, dx=1e-4)

    return optimize.root(lambda z: eq(z), [0.5])['x']


def func_opt(n_n, n_p, X):
    """
    X = [f, rho_0, rho_c, mu_c]
    """
    f, rho_0, rho_c = X
    n_in = np.array([f, n_n, n_p])
    mu_c = mu_ch(n_in, rho_0, rho_c)[0]
    eq1 = derivative(lambda z: eos.E_rho(np.array([z, n_n, n_p]), rho_0, rho_c, mu_c, C), f, dx=1e-4)
    eq2 = derivative(lambda z: eos.E_rho(n_in, z, rho_c, mu_c, C), rho_0, dx=1e-4, order=3)
    eq3 = derivative(lambda z: eos.E_rho(n_in, rho_0, z, mu_c, C), rho_c, dx=1e-4, order=3)
    return [eq1, eq2, eq3]

def sol(n_n, n_p, init=[0., 0., 0.]):
    return optimize.leastsq(lambda z: func_opt(n_n, n_p, z), init)[0]

nrange = np.linspace(10., 16.0, 100)
solutions = []
res = [0.0, 0., 0.]
for n in nrange:
    res = sol(n, 0., init=res)
    print(n, res)
    solutions.append(res)

lines = plt.plot(nrange, solutions)
plt.legend(lines, ['f', 'rho_0', 'rho_c'])
plt.show()
exit()
print(mu_ch(np.array([0.2, 0.6, 0.1]), 0., 0.,))

nrange = np.linspace(10., 16.0, 100)
print(rho_eq(8, 0.1))
plt.plot(nrange, list(map(lambda z: rho_eq(z, 0.), nrange)))
plt.show()
print(sqrt(C.Cr / C.eta_r(.2)) / C.M[0] / 2 * (0.1-0.5) / (C.m_rho * (1-0.2)))

exit()

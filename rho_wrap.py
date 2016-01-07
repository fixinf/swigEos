import eosWrap as eos
import numpy as np
from numpy import pi, sqrt
from scipy.misc import derivative
from scipy import optimize

def E_r(n_n, n_p, f, mu_c, C):
    return eos.E_rho(np.array([f, n_n, n_p]), mu_c, C)

def eq(n, C, z):
    n_p = z[0]
    f = z[1]
    mu_c = z[2]
    dx = 1e-4
    m_e = 0.5/135.
    m_mu = 105./135.
    n_in = np.array([f, n-n_p, n_p])
    mu_n = derivative(lambda z: E_r(z, n_p, f, mu_c, C), n-n_p, dx=dx)
    mu_p = derivative(lambda z: E_r(n - n_p, z, f, mu_c, C), n_p, dx=dx)
    mu_e = mu_n - mu_p
    n_rho = 2 * C.m_rho * C.M[0]**2 * sqrt(C.eta_r(f)) * C.phi_n(0, f)**2 / C.Cr * (1 - mu_c / (C.m_rho * C.phi_n(0,f)))
    n_rch = 0.
    if n_rho < n - 2 * n_p:
        r_c2 = (n - 2*n_p - n_rho) / (4 * C.m_rho * sqrt(C.eta_r(f)))
        n_rch += 2 * C.m_rho * C.phi_n(0, f) * r_c2
    n_e = 0.
    n_mu = 0.
    if mu_e > m_e:
        n_e += (mu_e**2 - m_e**2)**(3/2) / (3 * pi**2)
    if mu_e > m_mu:
        n_mu += (mu_e**2 - m_mu**2)**(3/2) / (3 * pi**2)

    eq_ch = n_p - n_e - n_mu - n_rch

    eq_f = derivative(lambda z: E_r(n - n_p, n_p, z, mu_e, C), f, dx=dx)

    return [[eq_ch, eq_f, mu_e - mu_c], [n_rho, n_rch]]

def eq_sol(n, C, init=[0., 0.2, 0.5]):
    res = optimize.leastsq(lambda z: eq(n, C, z)[0], init)[0]
    out = eq(n, C, res)[1]
    return [res, out]


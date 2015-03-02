import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join
from pylab import sqrt

C = Models.KVOR()
wr = Wrapper(C)
n = np.linspace(0., 8*wr.n0, 80, endpoint=False)
fs = []
mu = []
f = 0.
efRMF = []
for _n in n:
    n_f = np.array([_n, 0])
    f,= eos.f_eq(n_f, np.array([f]), 1, C)
    fs.append(f)
    mn = C.M[0]
    _mu = C.Co * _n / (mn**2 * C.eta_o(f)) + C.Cr * _n / (mn**2 * 4 * C.eta_r(f)) + sqrt(eos.p_f(_n)**2 + mn**2 * 
                                                      C.phi_n(0, f)**2)
    print _mu, eos.mu(np.array([f, _n, 0]), 2, C)
#     _mu = eos.mu(np.array([f, _n/2, _n/2]), 2, C)
    mu.append(_mu)
    
    efRMF.append(sqrt(eos.p_f(_n)**2 + mn**2 * C.phi_n(0, f)**2))
    
efRMF = np.array(efRMF)
mu = np.array(mu)
f0 = wr.f0_nm(n)
f1 = wr.f1_nm(n)
matsui_f1 = []
matsui_k = []
with open('Matsui_f1.csv', 'r') as f:
    for line in f:
        _k, _f1 = line.split(',')
        matsui_f1.append(float(_f1))
        matsui_k.append(float(_k))

                
print f1
print 1+f1/3

kf_plot = (0.16 / wr.n0)**(1./3) * np.array(map(eos.p_f, n))
# plt.plot(matsui_k, matsui_f1, kf_plot, f1)
# plt.show() 

ef = mu * (np.ones(f1.shape) + f1/3.)
kf = np.array(map(lambda z: eos.p_f(z), n))
meff = np.sqrt(ef**2 - kf**2)

plt.plot(kf_plot, 3 * (efRMF/mu - 1))
plt.plot(kf_plot, f1)
# plt.plot(matsui_k, matsui_f1)
plt.show()

plt.plot(n/wr.n0, meff/C.M[0])
# plt.plot(n/wr.n0, ef/C.M[0])
plt.plot(n/wr.n0, map(lambda z: C.phi_n(0, z), fs))
plt.show()

E = wr.Eneutr(n)
P = wr.P_N(n)/wr.const/wr.m_pi**4
vs = np.gradient(P)/np.gradient(E)
vs_lp = []
for i, _n in enumerate(n):
    vs_lp.append(kf[i]**2 / (3 * efRMF[i] * mu[i]) * (1 + f0[i]))

plt.plot(n/wr.n0, vs)
plt.plot(n/wr.n0, vs_lp)
plt.show()
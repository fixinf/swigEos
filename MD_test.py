import matplotlib
from scipy.misc.common import derivative
from time import sleep
from math import exp
matplotlib.use('QT4Agg')
import numpy as np
import matplotlib.pyplot as plt
import Models
from Wrapper import Wrapper
from scipy import optimize
import eosWrap as eos

C = eos.KVOR_MD()
wr = Wrapper(C)

KV = eos.KVOR()
wrKV = Wrapper(KV)

frange = np.linspace(0, 1, 1000)
n = np.linspace(0, 1*wr.n0, 100)

C.b_sigma = 0.23806
C.b_om = 0*0.953542
C.b_rho = 0.462841

# C.om_prime = 0.112492
C.om_prime = derivative(KV.eta_o, C.f0, dx=1e-3)
# print C.om_prime
# C.rho_prime = -0.391677
C.rho_prime = derivative(KV.eta_r, C.f0, dx=1e-3)

# plt.plot(frange, map(lambda z: derivative(C.eta_o, z, dx=1e-3), frange))
# plt.show()

print derivative(KV.eta_r, C.f0, dx=1e-3), derivative(KV.eta_o, C.f0, dx=1e-3)
# exit()

# fig, ax = plt.subplots(3, 1)



# ax[0].plot(frange, map(C.eta_r, frange), frange, map(KV.eta_r, frange))
# ax[0].set_ylim([0, 5])
# ax[1].plot(frange, map(C.eta_o, frange), frange, map(KV.eta_o, frange))
# ax[2].plot(frange, map(C.eta_s, frange), frange, map(KV.eta_s, frange))
# plt.show()

cAPR = 1.0/0.16
nAPR=cAPR*np.array([0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.32, 0.4, 0.48, 0.56, 0.64, 0.8, 0.96])
APR = np.array([0.0, -6.48, -12.13, -15.04, -16.00, -15.09, -12.88, -5.03, 2.13, 15.46, 
             34.39, 58.35, 121.25, 204.02])

APR_N = np.array([0., 4.45, 6.45 , 9.65, 13.29, 17.94, 22.92, 27.49, 38.82, 54.95, 75.13, 99.75, 127.58, 205.34, 305.87])
nAPR_N = cAPR* np.array([0, 0.02, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24,  0.32,  0.4, 0.48, 0.56, 0.64, 0.8, 0.96]) 

fig, ax = plt.subplots(2,2)
es, fs = wr.Esymm(n, ret_f=1)
esKV, fsKV = wrKV.Esymm(n, ret_f=1)
en, fn = wr.Eneutr(n, ret_f=1)
enKV, fnKV = wrKV.Eneutr(n, ret_f=1)


ax[0,0].plot(nAPR, APR, n/wr.n0, (es/n - C.M[0])*135, n/wr.n0, (esKV/n - C.M[0])*135)
for b in np.linspace(-1, 1, 20):
    C.b_rho = b
    en, fn = wr.Eneutr(n, ret_f=1)
    ax[1,0].plot(n/wr.n0, (en/n - C.M[0])*135)
    ax[1,1].plot(n/wr.n0, fn)

ax[1,0].plot(nAPR_N, APR_N, n/wr.n0, (enKV/n - C.M[0])*135)    
ax[0,0].set_xlim([0,1])
ax[1,0].set_xlim([0,1])
ax[0,0].set_ylim([-16, 2])
ax[1,0].set_ylim([0, 20])

ax[0,1].plot(n/wr.n0, fs, n/wr.n0, fsKV)
ax[1,1].plot(n/wr.n0, fn, n/wr.n0, fnKV)

plt.show()

wr.dumpVs()
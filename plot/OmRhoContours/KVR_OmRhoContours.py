#!/usr/bin/python
import eosWrap as eos
import matplotlib
from scipy.misc.common import derivative
matplotlib.use('QT4Agg')
import matplotlib.pyplot as plt
from Wrapper import Wrapper
import numpy as np
import os
from pylab import pause
from scipy import interpolate
imgfolder = os.path.join('..','img')

C = eos.KVOR_mod2()

C.SetHyperConstants(2)
wr = Wrapper(C)
n = np.linspace(0.0, 4.0, 2000)

uX = []
uY = []


C.Csp = 1

C.Cs = 179.56233875157545
C.Co =  87.59963973682763
C.Cr = 100.63642242792484
C.b = 0.007734608051455927
C.c = 0.0003446178665624873

C.rho_f = 0.48
C.rho_a = 0*1000.
C.omega_f = 0.48
C.omega_a = 0*100.
C.omega_kind = 1
C.rho_power = 1.
f0 = 0.2
C.rho_kind = 1

suffix = 'KVR_aom=%.2f'%C.omega_a

C.beta = 2.9
C.gamma = 2.
C.alpha = 0.0
C.SetHyperConstants(2)

npoints = 1000

with open('../klahnUpper2', 'r') as f:
    for line in f:
        x,y = line.split()
        uX.append(float(x)/0.16)
        uY.append(float(y))

# plt.semilogy(uX, uY)
x = 1.00
# print (wr.Esymm([x*wr.n0])/(x*wr.n0) - C.M[0])*wr.m_pi
# print eos.K(wr.n0, C)
# print eos.J(wr.n0, C)
# print eos.f_eq(np.array([wr.n0/2, wr.n0/2]), np.array([0.0]), 1, C)
fig, ax = plt.subplots()

UKlahnY = []
UKlahnX = []
with open('../klahnUpper', 'r') as f:
    for line in f:
        _n, p = line.strip().split()
        UKlahnY.append(float(p))
        UKlahnX.append(float(_n)/0.16)
ax.semilogy(UKlahnX, UKlahnY, c = 'grey')

LKlahnX = []
LKlahnY = []
with open('../klahnLower', 'r') as f:
    for line in f:
        _n, p = line.strip().split()
        LKlahnY.append(float(p))
        LKlahnX.append(float(_n)/0.16)
ax.semilogy(LKlahnX, LKlahnY, c = 'grey')

print C.f0
# wr.solve()
C.f0 = f0
wr.solve(f0 = C.f0, E0=-15.8, K0 = 250, J0=28.0, iter = 3000)
print C.f0
# exit()
# set(params)
C.SetHyperConstants(2)
# print C.alpha, C.d

flist = [f0]#[0.195, 0.26, 0.27, 0.28]
lines = []
labels = []
n_star = []
mlist = []
Vs = []
rho_a_list = [10**i for i in np.linspace(0, 4., 20)]
omega_a_list = [10**i for i in np.linspace(0, 4., 20)]
fContour = open('KVRContour.dat', 'w')
for rho_a in rho_a_list:
    for omega_a in omega_a_list:
        global n_star
    #     for _f in np.linspace(C.f0, f):
        print C.Cs, C.Co, C.Cr, C.b, C.c    
        print 'J = ', eos.K(wr.n0, C)
        print 'K = ',eos.J(wr.n0, C)
        print 'K\' = ',-3*wr.n0*(derivative(lambda z: eos.K(z, C), wr.n0, dx=1e-3, order=3) - 
                       2*eos.K(wr.n0,C)/wr.n0)
        print 'L = ', 3*wr.n0*derivative(lambda z: eos.J(z, C), wr.n0, dx=1e-3)
    #     pause()
        C.omega_a = omega_a
        C.rho_a = rho_a
        wr.reset(npoints=npoints, iter=100)
        wr.setDriver()
        n_star, M, R = wr.stars(npoints=200)
        mlist.append(M)
        fContour.write('%f  '%(np.max(M)))
    fContour.write("\n")
f.close()

# with open('masses_noscaling.dat', 'w') as f:
    # for i, _n in enumerate(n_star):
    #     f.write('%f   '%(_n/wr.n0))
    #     for _m in mlist[:, i]:
    #         print _m, mlast
    #         if mlast > _m:
    #             print 'breaking!'
    #             break
    #         f.write('%f   '%_m)
    #         mlast = _m
    #     f.write('\n')
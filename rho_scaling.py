#!/usr/bin/python
import eosWrap as eos
import matplotlib
from scipy.misc.common import derivative
from scipy.constants.constants import pi
from numpy import sqrt
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

# C.Cs = 179.56233875157545
# C.Co =  87.59963973682763
# C.Cr = 100.63642242792484
# C.b = 0.007734608051455927
# C.c = 0.0003446178665624873

C.Cs = 227.825989
C.Co =  134.882595
C.Cr = 93.199075
C.b = 0.005592
C.c = -0.012123

# C.alpha = 0.85
# C.z = 0.65
#     
# C.omega_a = 6.45
# C.omega_f = 0.53
# 
# C.phi_a = -0.85
# C.phi_f = 0.28
#     
# C.d = -5.5
#     
# C.phi_gamma = 3.0
# C.phi_z = 3.5
#  
# C.sprime = 0
# C.Csp = 380.0
#     
# C.beta = 0.8
# C.gamma = 9.0
# C.exp_alpha = 0.05
# C.omega_c = -15000.0

C.rho_f = 0.4
C.rho_a = 100
# print C.rho_a
C.beta = 1.0
C.gamma = 5.5
C.f0 = 0.195

#C.rho_kind = 0

def eval_rho_a(f_max, f_rho, f0, beta, gamma):
    return ((1 + beta*f_max**2)/(1 + beta*f0**2))**gamma / (f_max-f_rho)**3

# C.rho_a =  -eval_rho_a(1.1, C.rho_f, C.f0, C.beta, C.gamma)
print C.rho_a 
pause(1)
C.SetHyperConstants(2)
print [i for i in C.X_s]
C.set_xs(np.array([0.0, 0.0, -30.0, 30.0, 30., 30., -18., -18.]))
print [i for i in C.X_s]

wr.solve(f0=C.f0,iter=3000)
C.SetHyperConstants(2)
C.phi_meson = 1
hyper=0
npoints=1000
wr.reset(hyper=hyper, npoints=npoints)

rho = []
for r in wr.rho:
    sum = 0
    for _n in r[1:]:
        sum += _n
    rho.append(r[1:]/sum)


def diffJ(n, f=None, f_init=0):
    if f is None:
        f = eos.f_eq(np.array([n/2, n/2]), np.array([f_init]), 1, C)[0]
#     print f_eq
    ES = eos.EBind(np.array([f, n/2, n/2]), C)
    f = eos.f_eq(np.array([n, 0]), np.array([f]), 1, C)[0]
    EN = eos.EBind(np.array([f, n, 0]), C)
    return EN - ES

def anJ(n, f=None, f_init=0):
    if f is None:
        f = eos.f_eq(np.array([n/2, n/2]), np.array([f_init]), 1, C)[0]
    res = C.Cr * n / (8* C.M[0]**2 * C.eta_r(f))
    p_f = (3 * pi **2 * n/2)**(1.0/3.0)
    m_eff = C.M[0]*C.phi_n(0, f)
    res += pi**2 * n / (4*p_f * sqrt(p_f**2 + m_eff**2))
#     print 'res=', res
    return res*wr.m_pi

dJ = []
aJ = []
nJ = []
eta_r=[]
f_s_list=[]
f_eq = 0.0
for n in wr.n:
    global f_eq
    f_eq = eos.f_eq(np.array([n/2, n/2]), np.array([f_eq]), 1, C)[0]
    dJ.append(diffJ(n, f=f_eq))
    aJ.append(anJ(n, f=f_eq))
    nJ.append(eos.J(n, C, f_eq))
    eta_r.append(C.eta_r(f_eq))
    f_s_list.append(f_eq)

# #------------ Plot three ways ----------------
# fig= plt.figure(figsize=(6,6))
# Jlist = np.array([dJ, aJ, nJ]).transpose()
# lines = plt.plot(wr.n/wr.n0, Jlist)
# plt.ylim([0.0, 300.0])
# plt.legend(lines, ['difference','analytic', 'numeric'], loc=0)
# plt.xlabel(r'$n/n_0$')
# plt.ylabel(r'$J(n) [MeV]$')
# plt.show()
# #-------------- end -------------------      

# #-------------- Plot the pole ---------------#
# fig= plt.figure(figsize=(6,6))
# Jlist = np.array([aJ, nJ]).transpose()
# lines = plt.plot(wr.n/wr.n0, Jlist)
# plt.ylim([0.0, 300.0])
# plt.legend(lines, ['analytic', 'numeric'], loc=0)
# plt.xlabel(r'$n/n_0$', fontsize=18)
# plt.ylabel(r'$J(n) [MeV]$',fontsize=18)
# plt.show()
# #-------------- end -------------------#

# #-------------- Plot \eta_\rho ---------------#
# fig= plt.figure(figsize=(6,6))
# Jlist = np.array([aJ, nJ]).transpose()
# lines = plt.plot(wr.n/wr.n0, eta_r, 
#                  wr.n/wr.n0, map(C.eta_r, wr.rho[:,0]))
# plt.ylim([-1.0, 5.0])
# plt.legend(lines, ['symmetric', 'NS matter'], loc=0)
# plt.xlabel(r'$n/n_0$',fontsize=18)
# plt.ylabel(r'$\eta_\rho(f_{eq}(n))$',fontsize=18)
# plt.show()
# #-------------- end -------------------#

# #-------------- Plot f_eq ---------------#
# fig= plt.figure(figsize=(6,6))
# Jlist = np.array([aJ, nJ]).transpose()
# lines = plt.plot(wr.n/wr.n0, f_s_list, 
#                  wr.n/wr.n0, wr.rho[:,0])
# plt.ylim([0.0, 1.0])
# plt.legend(lines, ['symmetric', 'NS matter'], loc=0)
# plt.xlabel(r'$n/n_0$', fontsize=18)
# plt.ylabel(r'$f_{eq}$', fontsize=18)
# plt.show()
# #-------------- end -------------------#

# #-------------- Plot J(n) for two c_omega ---------------#
# c_list = [-1000.0, -350000.0]
# fig= plt.figure(figsize=(6,6))
# lines = []
# lines2 = []
# lslist = ['--', '-']
# for j, c in enumerate(c_list):
#     dJ = []
#     aJ = []
#     nJ = []
#     eta_r=[]
#     f_s_list=[]
#     f_eq = 0.0
#     C.omega_c = c
#     for n in wr.n:
#         global f_eq
#         f_eq = eos.f_eq(np.array([n/2, n/2]), np.array([f_eq]), 1, C)[0]
#         dJ.append(diffJ(n, f=f_eq))
#         aJ.append(anJ(n, f=f_eq))
#         nJ.append(eos.J(n, C, f_eq))
#         eta_r.append(C.eta_r(f_eq))
#         f_s_list.append(f_eq)
#     line, = plt.plot(wr.n/wr.n0, map(C.eta_r, f_s_list), c='g', ls=lslist[j])
#     wr.reset(npoints=1000)
#     line2, = plt.plot(wr.n/wr.n0, map(C.eta_r, wr.rho[:, 0]), c='b', ls=lslist[j]) 
#     lines.append(line)
#     lines2.append(line2)
# 
# plt.ylim([-1.0, 5.0])
# print lines+lines2
# plt.legend(lines+lines2, 
#            ['old symm', 'new symm','old NS', 'new NS'], loc=0)
# # plt.legend(lines2, ['new symm', 'new NS'], loc=0)
# plt.xlabel(r'$n/n_0$', fontsize=18)
# plt.ylabel(r'$\eta_\rho(f_{eq}(n))$', fontsize=18)
# plt.show()
# #-------------- end -------------------#

# #-------------- Plot hyper for two a_rho ---------------#
# colors = ['b', 'r', 'g', 'y', 'cyan', 'm', 'k', 'black']
# ar_list = [0.0, -70.466]
# fig= plt.figure(figsize=(6,6))
# lines = []
# lines2 = []
# mlist = []
# lslist = ['--', '-']
# for j, a in enumerate(ar_list):
#     rho = []
#     C.rho_a = a
#     plt.gca().set_color_cycle(colors)
#     wr.reset(hyper=1, npoints=npoints)
#     for r in wr.rho:
#         sum = 0
#         for _n in r[1:]:
#             sum += _n
#         rho.append(r[1:]/sum)
#     lines=plt.plot(wr.n/wr.n0, rho, ls=lslist[j])
#     wr.setDriver()
#     N, M, R = wr.stars(npoints=200)
#     mlist.append(np.max(M))
# #     line2, = plt.plot(wr.n/wr.n0, map(C.eta_r, wr.rho[:, 0]), c='b', ls=lslist[j]) 
# #     lines.append(line)
# #     lines2.append(line2)
# 
# plt.ylim([0.0, 1.0])
# # plt.legend(lines,['old', 'new'], loc=0)
# # plt.legend(lines2, ['new symm', 'new NS'], loc=0)
# plt.legend(lines, ['n','p',r'$\Lambda$',r'$\Sigma^-$',r'$\Sigma^0$',
#            r'$\Sigma^+$', r'$\Xi^-$', r'$\Xi^0$'],
#            loc=0)
# plt.xlabel(r'$n/n_0$', fontsize=18)
# plt.ylabel(r'$n_i$', fontsize=18)
# print mlist
# plt.show()
# #-------------- end -------------------#

fig, ax = plt.subplots(2,4)
ax[0,0].plot(wr.n/wr.n0, rho, wr.n/wr.n0, [0.14 for i in wr.n])
E, P, N = wr.EPN()
ax[0,1].plot(N[1:]/wr.n0, np.diff(P)/np.diff(E))
ax[0,2].semilogy(wr.n[1:]/wr.n0, wr.Psymm(wr.n))
ax[0,3].plot(wr.n/ wr.n0, dJ)
ax[0,3].plot(wr.n/wr.n0, nJ)
ax[0,3].plot(wr.n/wr.n0, aJ)
ax[0,3].set_ylim([0.0, 300.0])
ax[1,0].plot(wr.n/wr.n0, map(lambda z: C.eta_r(z), wr.rho[:,0]))
ax[1,0].plot(wr.n/wr.n0, eta_r)
ax[1,1].plot(wr.n/wr.n0, f_s_list)
ax[1,1].plot(wr.n/ wr.n0, wr.rho[:, 0])
wr.setDriver()
N, M, R = wr.stars(npoints=200)
print max(M)
ax[1,2].plot(N/wr.n0, M)
plt.show()
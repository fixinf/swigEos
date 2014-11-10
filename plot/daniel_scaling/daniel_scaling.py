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
import Models
imgfolder = os.path.join('..','img')

C = Models.myMod()

C.SetHyperConstants(2)
wr = Wrapper(C)
n = np.linspace(0.0, 4.0, 2000)

uX = []
uY = []

f0 = 0.27

C.Csp = 1


suffix = 'Mod/KVOR_rhoa_%.2f'%C.rho_a

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
wr.solve(f0 = C.f0, iter = 3000)
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
for f in flist:
    global n_star
#     for _f in np.linspace(C.f0, f):
    wr.solve(f0=f, iter=3000)
    print C.Cs, C.Co, C.Cr, C.b, C.c    
    print 'J = ', eos.K(wr.n0, C)
    print 'K = ',eos.J(wr.n0, C)
    print 'K\' = ',-3*wr.n0*(derivative(lambda z: eos.K(z, C), wr.n0, dx=1e-3, order=3) - 
                   2*eos.K(wr.n0,C)/wr.n0)
    print 'L = ', 3*wr.n0*derivative(lambda z: eos.J(z, C), wr.n0, dx=1e-3)
#     pause()
    C.omega_c = -1000.0
    l1, = ax.plot(n[:]/wr.n0, wr.Psymm(n), c='r', ls='--')
    C.omega_c = -15000.0
    l2, = ax.plot(n[:]/wr.n0, wr.Psymm(n), c='r')
    plt.xlabel(r'$n/n_0$', fontsize=18)
    plt.ylabel(r'$P_{symm}[MeV]$', fontsize=18)
    plt.ylim([1.0, 1000.0])
    plt.legend([l1, l2], ['old', 'new'], loc=0)
    plt.show()
    lines.append(wr.Psymm(n))
    labels.append(str(C.phi_n(0, f)))
    wr.reset(npoints=npoints, iter=100)
    
    wr.setDriver()
    n_star, M, R = wr.stars(npoints=200)
    mlist.append(M)

_P = wr.Psymm(n)
_E = wr.Esymm(n)*wr.m_pi**4*wr.const #conversion to MeV/fm^3
dE = np.diff(_E)
dP = np.diff(_P)

wr.reset(npoints=1000, nmax=25.)
_ENS, _PNS, _N = wr.EPN()
plt.plot(_N/wr.n0, _PNS)
plt.show()
# plt.plot(_N/wr.n0, _PNS)
# plt.show()
Vs = dP/dE[1:]
VsNS = np.diff(_PNS)/(np.diff(_ENS))
plt.plot(_N[1:]/wr.n0, VsNS)
plt.plot(n[2:]/wr.n0, Vs)
plt.show()

with open(suffix+'vs_scaling%.2f.dat'%flist[0], 'w') as f:
    for i, _n in enumerate(n[2:]):
        f.write('%f    %f \n' %(_n/wr.n0, Vs[i]))
        
with open (suffix+'vsNS_scaling%.2f.dat'%flist[0], 'w') as f:
    for i, _n in enumerate(_N[1:]):
        f.write('%f    %f \n' %(_n/wr.n0, VsNS[i]))

lines = np.array(lines)
mlist = np.array(mlist)

with open(suffix+'pressure_scaling%.2f.dat'%flist[0], 'w') as f:
    for i, _n in enumerate(n[1:]):
        f.write('%f   '%(_n/wr.n0))
        for _p in lines[:,i]:
            f.write('%f   '%_p)
        f.write('\n')


for j,m in enumerate(mlist):
    with open(suffix+'mass_scaling_%.2f.dat'%flist[j],'w') as f:
        for i, _n in enumerate(n_star):
            print m[i]
            print max(m)
            if i > np.argmax(m):
                break
            f.write('%f  %f\n'%(_n/wr.n0, m[i]))
    MofN = interpolate.interp1d(n_star, m)



rho = []
sum = []
for r in wr.rho[:,1:]:
    rho.append(r/np.sum(r))
    
with open(suffix+'n_scaling_%.2f.dat'%flist[j], 'w') as f:
    for i, r in enumerate(rho):
        f.write("%f  " % (wr.n[i]/wr.n0))
        for k in r:
            f.write("%f   " % k)
        f.write("\n")
        
with open(suffix+'f_%.2f.dat'%flist[j], 'w') as f:
    for i, _n in enumerate(wr.n):
        f.write('%f %f\n'% (_n/wr.n0, wr.rho[i][0]))

with open(suffix+'eta_rho_n_%.2f.dat'%flist[j], 'w') as f:
    for i, _n in enumerate(wr.n):
        f.write('%f %f\n'% (_n/wr.n0, C.eta_r(wr.rho[i][0])))
        
with open(suffix+'eta_rho_f_%.2f.dat'%flist[j], 'w') as f:
    for _f in np.linspace(0.0, 1.0, 100):
        f.write('%f %f \n' % (_f, C.eta_r(_f)))

plt.plot(wr.n/wr.n0, rho, wr.n/wr.n0, [0.14 for i in wr.n])
plt.show()
rho = np.array(rho)
i_DU = np.argmin(np.abs(rho[:,1] - np.array([0.14 for r in rho])))
n_DU = wr.n[i_DU]
M_DU = MofN(n_DU)

print n_DU/wr.n0, M_DU



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
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
from pylab import log, sqrt, pi
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
C.rho_a = 0*100.
C.omega_f = 0.2
C.omega_kind = 1
C.omega_a = 1000.
C.rho_power = 2.
f0 = 0.195
C.rho_kind = 1

suffix = 'KVOR_fom=%.2f_aom=%.2f'%(C.omega_f, C.omega_a)

C.beta =  1.0
C.gamma = 11.
C.alpha = 1.0
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
def KP(n, m):
    p_f = eos.p_f(n)
    res = -(3./8.)*m**4 * log(sqrt(m**2))
    res += (3./8.)*m**4*log(p_f+sqrt(m**2+p_f**2))
    res += sqrt(m**2+p_f**2)* (-((3*m**2 *p_f)/8.)+p_f**(3)/4.0)
    return res/(3*pi**2)

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
    l1, = ax.plot(n[1:]/wr.n0, wr.Psymm(n), c='r', ls='--')
    l2, = ax.plot(n[1:]/wr.n0, wr.Psymm(n), c='r')
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



terms = []
f_eq_n = 0.
def E_N(n):
    global f_eq_n
    f_eq, = eos.f_eq(np.array([n, 0.]), np.array([f_eq_n]), 1, C)
    f_eq_n = f_eq
    return eos._E(np.array([f_eq, n, 0.]), C)

EN = map(E_N, wr.n)
dn = wr.n[1]-wr.n[0]
PN = wr.n[1:]*np.diff(EN)/dn - EN[1:]

f_eq = 0.
with open(suffix+"terms_%.2f.dat"%flist[0], 'w') as f:
    for i, _n in enumerate(wr.n[:]):
        n_n = _n
        n_p = 0.
        f_eq, = eos.f_eq(np.array([n_n, 0.]), np.array([f_eq]), 1, C)
        omega = C.Co * _n**2/(2*C.M[0]**2 * C.eta_o(f_eq))
        sigma = C.M[0]**4 * f_eq**2 / (2 * C.Cs) * C.eta_s(f_eq)
        _rho = C.Cr * (n_n - n_p)**2 / (8 * C.M[0]**2 * C.eta_r(f_eq))
        K = eos.kineticInt(n_n, C.M[0]*C.phi_n(0, f_eq), f_eq)
        dn = 1e-3
#         f_eq, = eos.f_eq(np.array([n_n + dn, 0.]), np.array([f_eq]), 1, C)
        dK = eos.kineticInt(n_n + dn, C.M[0]*C.phi_n(0, f_eq), f_eq)
#         f_eq, = eos.f_eq(np.array([n_n-dn, 0.]), np.array([f_eq]), 1, C)
        dK -= eos.kineticInt(n_n -dn, C.M[0]*C.phi_n(0, f_eq), f_eq)
        dK /= dn*2
        K += eos.kineticInt(n_p, C.M[0]*C.phi_n(0, f_eq), f_eq)
        K1 = KP(n_n, M[0]*C.phi_n(0, f_eq))
#         K1 += KP(n_p, M[0]*C.phi_n(0, f_eq))
        print n_n*dK - K, K1
#         print f_eq, n_n, n_p, PN[i], (n_n*dK - K + omega - sigma + _rho), EN[i], (K + omega + sigma + _rho)
        f.write('%f %f %f %f %f %f \n'%(_n/wr.n0, sigma, omega, _rho, K, n_n*dK - K))
        terms.append([sigma, omega, _rho, K, n_n*dK - K])

linest = plt.plot(wr.n/wr.n0, terms)
plt.legend(linest, ['sigma', 'omega', 'rho', 'K', 'K_2'], loc=0)
plt.show()

_P = wr.Psymm(n)
_E = wr.Esymm(n)*wr.m_pi**4*wr.const #conversion to MeV/fm^3
dE = np.diff(_E)
dP = np.diff(_P)

wr.reset(npoints=400)
_ENS, _PNS, _N = wr.EPN()

with open(suffix+"PE_%.2f.dat"%flist[0], 'w') as f:
    for i,p in enumerate(_PNS):
        f.write('%f %f \n'%(_ENS[i], _PNS[i]))

plt.plot(_ENS, _PNS)
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
            print _n, m[i]
            print max(m)
            if i > np.argmax(m):
                print _n
                print wr.n[np.argmin(abs(wr.n - [_n for i in wr.n]))]
                print '%e'%(wr.E[np.argmin(abs(wr.n - [_n for i in wr.n]))]*wr.m_pi**4)
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

with open(suffix+'eta_om_n_%.2f.dat'%flist[j], 'w') as f:
    for i, _n in enumerate(wr.n):
        f.write('%f %f\n'% (_n/wr.n0, C.eta_o(wr.rho[i][0])))

with open(suffix+'eta_om_f_%.2f.dat'%flist[j], 'w') as f:
    for _f in np.linspace(0.0, 1.0, 100):
        f.write('%f %f \n' % (_f, C.eta_o(_f)))

plt.plot(wr.n/wr.n0, rho, wr.n/wr.n0, [0.14 for i in wr.n])
plt.show()
rho = np.array(rho)
i_DU = np.argmin(np.abs(rho[:,1] - np.array([0.14 for r in rho])))
n_DU = wr.n[i_DU]
M_DU = MofN(n_DU)

nmax = np.argmax(m)

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
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
C1 = Models.myMod()
C2 = Models.myMod()
wr1 = Wrapper(C1)
wr2 = Wrapper(C2)

frange = np.linspace(0., 1., 1000)


# C1.rho_f = 1.
# C1.gamma = 0.

print C2.eta_r(0.59)
print derivative(lambda z: C1.eta_r(z), 0.27, dx=1e-3)
# exit()

# TEST OF RHO KIND 6
C2.rho_kind = 6
C2.rho_a = 0.6
C2.rho_sat_a = 20
C2.rho_sat_f1 = 0.35
C2.beta = 3.4
C2.gamma = 0
C2.rho_sat_val = 9.99673914
C2.rho_f = 0.65309008
# C2.rho_f_low = 0.29
# C2.rho_a_low = 25.
C2.omega_a_low = 0
def f(x):
    C2.rho_f = x[0]
    C2.rho_sat_val = x[1]
#     return [C2.eta_r(0.27) - 1,
#             derivative(C2.eta_r, 0.27, dx=1e-3) - 
#             derivative(C1.eta_r, 0.27, dx=1e-3)]

    return [C2.eta_r(0.27)-1, derivative(C2.eta_r, 0.27, dx=1e-3) - 
            derivative(C1.eta_r, 0.27, dx=1e-3)]

print optimize.broyden1(f, [C2.rho_f, C2.rho_sat_val])

print [derivative(C2.eta_r, 0.27, dx=1e-3) - 
            derivative(C1.eta_r, 0.27, dx=1e-3)]

etar_old = map(C2.eta_r, frange)

C2.drho = derivative(C2.eta_r, C2.f0, dx=1e-3)
C2.d2rho = derivative(C2.eta_r, C2.f0, dx=1e-3, n=2)


C2.rho_a_low = 1.1

C2.rho_b_low = -0.5
# C2.rho_b_low = -C2.f0/5.*(C2.drho - C2.d2rho*C2.f0 + 3*C2.rho_a_low/C2.f0)

C2.rho_kind = 7


# wr1=wr2

npoints = 100

S = np.array(map(C2.s, frange))
# plt.plot(frange, S, frange, frange*(1-S) + S*C2.rho_sat_f2)
# plt.show()

# plt.plot(frange, map(C1.eta_o, frange))
# plt.plot(frange, map(C2.eta_o, frange))
# plt.show()

plt.plot(frange, map(C1.eta_r, frange))
plt.plot(frange, map(C2.eta_r, frange))

plt.plot(frange, etar_old)
plt.show()

order = 3
plt.plot(frange, map(lambda z: derivative(C2.eta_r, z, dx=1e-3, n=order, order = 5), frange))
plt.plot(frange, map(lambda z: derivative(C2.eta_o, z, dx=1e-3, n=order, order=5), frange))
plt.plot(frange, map(lambda z: derivative(C2.eta_s, z, dx=1e-3, n=order, order=5), frange))
plt.show()


n = np.linspace(0, 1*wr1.n0, 100)
# C2.Cr = 91
# for J in np.linspace(30., 32., 2):
wr2.solve(f0=C2.f0, K0=240., J0=30)
print wr2.J(), wr2.L()

data = np.loadtxt('/home/const/GrabbedFigures/HebelerSchwenk2014/HebelerSchwenk_Out.dat',
                   skiprows=1)

fig, ax = plt.subplots(2,1)
good = []
for b in np.linspace(0, 2.5, 10):
    flag = 0
    C2.rho_a_low = b
#     C2.rho_b_low = -C2.f0/5.*(C2.drho - C2.d2rho*C2.f0 + 3*C2.rho_a_low/C2.f0)
    wr2.reset(hyper=0, nmin=0., nmax=1*wr2.n0, npoints=100)
    vsNs = np.diff(wr2.P)/np.diff(wr2.E)
    vlast = 0
    for v in vsNs:
        if v < vlast:
            flag = 1
            break
        vlast = v
    if not flag:
        good.append(C2.rho_a_low) 
        ax[0].plot(n/wr1.n0, 135*(wr2.Eneutr(n)/n - C1.M[0]))
        ax[1].plot(wr2.n[1:]/wr2.n0, vsNs)
ax[0].plot(n/wr1.n0, 135*(wr1.Eneutr(n)/n - C1.M[0]), c='red',lw=2)
ax[0].plot(data[:,0], data[:, 1:],c='black',lw='3')
ax[0].set_xlim([0,1])
print good
plt.show()
# exit()
wr1=wr2
# wr1.C.phi_meson = 1
wr1.dumpVs()
wr1.testDU()
wr1.setDriver()
n, m, r, mg = wr1.stars()
print 'M_NS_MAX =', max(m)
sleep(3)
fNS = wr1.rho[:,0]
contribNS = wr1.getContrib()
nNS = wr1.n
print max(fNS)
# exit()
wr1.C.phi_meson=0
wr1.reset(npoints=npoints, timeout=2, hyper=1)
wr1.setDriver()
n, m, r, mg = wr1.stars()
print max(m)
# exit()
fig, ax = plt.subplots(2, 2)
ax[0,0].plot(wr1.n/wr1.n0, wr1.getContrib(), nNS/wr1.n0, contribNS)
ax[1,0].plot(wr1.n/wr1.n0, wr1.concentrations(), wr1.n/wr1.n0, wr1.rho[:,0], nNS/wr1.n0, fNS)
ax[0,1].plot(wr1.rho[:,0], map(wr1.C.eta_r, wr1.rho[:,0]))

plt.show()



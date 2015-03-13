import matplotlib
from scipy.misc.common import derivative
from time import sleep
from math import exp
from cmath import sqrt
matplotlib.use('QT4Agg')
from matplotlib.widgets import Slider
import numpy as np
import matplotlib.pyplot as plt
import Models
from Wrapper import Wrapper
from scipy import optimize
C1 = Models.myMod()
C2 = Models.myMod()
wr1 = Wrapper(C1)
wr2 = Wrapper(C2)

# C1.rho_f = 1.
# C1.gamma = 0.

print C2.eta_r(0.59)
print derivative(lambda z: C1.eta_r(z), 0.27, dx=1e-3)
# exit()

rho_f = C2.rho_f
def f(x):
    C2.rho_a = x[0]
    C2.rho_sat_val = x[1]
    C2.beta = x[2]
    C2.rho_a0 = x[3]
    res =[C2.eta_r(0.27)-1, derivative(C2.eta_r, 0.27, dx=1e-3) - 
            derivative(C1.eta_r, 0.27, dx=1e-3), C2.eta_r(C2.rho_f) - C1.eta_r(rho_f), C2.eta_r(1) - 1.5]
    res = np.array(res)
    return np.sum(res*res)

npoints = 100

frange = np.linspace(0., 1, 1000)
lines = []
labels = []
dlist = np.linspace(0, C2.rho_a0, 10)


for d in dlist:
    C2.rho_a0 = d
    res = optimize.minimize(f, [C2.rho_f, C2.rho_sat_val, C2.beta, C2.rho_a0]).x

fig, ax = plt.subplots()

ax.plot(frange, map(C2.eta_r, frange))

# wr2.testPodsiedlowski()

# plt.show()
n_full = 100

wr2.reset(npoints=n_full, iter=200)
print res
# En1 = wr2.Eneutr(n_full)
# Pn1 = wr2.P_N(n_full)/wr2.const/wr2.m_pi**4
vs1 = np.diff(wr2.P)/np.diff(wr2.E)
np1 = wr2.concentrations()[:,1]
fns1 = wr2.rho[:,0]
etar1 = map(C2.eta_r, fns1)
 
# 
# C2.rho = 0
# C2.drho = 0
# C2.d2rho = 0

C2.rho_ld = 0

f1 = 0.4

C2.rho = C1.eta_r(f1)
C2.drho = derivative(C1.eta_r, f1, dx=1e-3)
C2.d2rho = derivative(C1.eta_r, f1, dx=1e-3, n=2)


print C2.rho, C2.drho, C2.d2rho


C2.rho_a_low = 0.04
C2.rho_b_low = -0.03

fsmall = np.linspace(0., 0.27, 20)

def fun_low(x):
    C2.rho_a_low = x[0]
    C2.rho_b_low = x[1]
    res = np.array(map(lambda z: C2.eta_r(z) - C1.eta_r(z), fsmall))
    return np.sum(res**2)

C2.rho_sat_f2 = f1

# print optimize.minimize(fun_low, [0,0])

C2.rho_tan_c = 0

print wr2.L(), wr2.J()



line, = ax.plot(frange, map(C2.eta_r, frange))
lines.append(line)
# labels.append('%.2f'%d)
# ax.legend(lines, labels, loc=0)
ax.plot(frange, map(C1.eta_r, frange))
ax.set_ylim([0,5])
ax.set_xlabel(r'$f$', fontsize=18)
ax.set_ylabel(r'$\eta_\rho(f)$', fontsize=18)
#   

# ax_a = plt.axes([0.25, 0.1, 0.65, 0.03])
# sl_a = Slider(ax_a, 'A', -20, 20, valinit=C2.rho_tan_a)
#   
# ax_b = plt.axes([0.25, 0.15, 0.65, 0.03])
# sl_b = Slider(ax_b, 'B', -20, 20, valinit=C2.rho_tan_b)
#   
# ax_c = plt.axes([0.25, 0.2, 0.65, 0.03])
# sl_c = Slider(ax_c, 'C', -20, 20, valinit=C2.rho_tan_c)
# 
# ax_al = plt.axes([0.25, 0.25, 0.65, 0.03])
# sl_al = Slider(ax_al, 'AL', -0.5, 0.5, valinit=0)
# 
# ax_bl = plt.axes([0.25, 0.3, 0.65, 0.03])
# sl_bl = Slider(ax_bl, 'BL', -0.5, 0.5, valinit=0)
# 
#  
# 
# def update(val):
#     C2.rho_tan_a = sl_a.val
#     C2.rho_tan_b = sl_b.val
#     C2.rho_tan_c = sl_c.val
#     C2.rho_a_low = sl_al.val
#     C2.rho_b_low = sl_bl.val
#     line.set_ydata(map(C2.eta_r, frange))
#     
# sl_a.on_changed(update)
# sl_b.on_changed(update)
# sl_c.on_changed(update)
# sl_al.on_changed(update)
# sl_bl.on_changed(update)
  
plt.show()

C3 = Models.KVOR()

l1, = plt.plot(frange, map(lambda z: (1-z)**2/sqrt(C2.eta_r(z)), frange))
lKV, = plt.plot(frange, map(lambda z: (1-z)**2/sqrt(C3.eta_r(z)), frange))
plt.legend([l1, lKV], ['Current scaling', 'KVOR'], loc=0)
plt.xlabel(r'$f$', fontsize=18)
plt.ylabel(r'$\chi_\rho(f)$', fontsize=18)

plt.show()

plt.plot(wr2.n/wr1.n0, wr1.ESbind(wr2.n), wr2.n/wr1.n0, wr1.ESbind(wr2.n))
plt.show()

wr2.solve(f0=C2.f0, K0=240., J0=30.)

wr2.reset(npoints=n_full)
# En2 = wr2.Eneutr(n_full)
# Pn2 = wr2.P_N(n_full)/wr2.const/wr2.m_pi**4
vs2 = np.diff(wr2.P)/np.diff(wr2.E)
np2 = wr2.concentrations()[:,1]
fns2 = wr2.rho[:,0]
etar2 = map(C2.eta_r, fns2)

wr2.solve(f0=C2.f0, K0=240., J0=30.)
print wr2.J(), wr2.L()

fig, ax = plt.subplots(2,2)
ax[0,0].plot(wr2.n[1:]/wr2.n0, vs1, wr2.n[1:]/wr2.n0, vs2)
ax[1,0].plot(wr2.n/wr2.n0, np1, wr2.n/wr2.n0, np2)
ax[0,1].plot(wr2.n/wr2.n0, fns1, wr2.n/wr2.n0, fns2)
ax[1,1].plot(wr2.n/wr2.n0, etar1, wr2.n/wr2.n0, etar2)
print C2.rho_tan_a, C2.rho_tan_b, C2.rho_tan_c
plt.show()

Es, fs = wr2.Esymm(wr2.n,ret_f=1)
plt.plot(wr2.n/wr2.n0, map(lambda z: 1./C2.eta_r(z),fs))
plt.plot(wr2.n/wr2.n0, map(lambda z: 1./C1.eta_r(z),fs))
plt.show()

n, J2 = wr2.dumpJ()
n, J1 = wr1.dumpJ()
plt.plot(n/wr1.n0, J1, n/wr1.n0, J2)
plt.show()

# wr2.reset()
# wr1.reset()
# lines = plt.plot(wr1.n/wr1.n0, wr1.rho[:,0], wr2.n/wr2.n0, wr2.rho[:,0])
# plt.legend(lines, ['old', 'new'], loc=0)
# plt.show()

wr1=wr2
# wr1.C.phi_meson = 1
# wr1.dumpVs()


data = np.loadtxt('/home/const/GrabbedFigures/HebelerSchwenk2014/HebelerSchwenk_Out.dat',
                   skiprows=1)
n = np.linspace(0, wr1.n0, 100)
# plt.plot(data[:,0], data[:, 1:],c='black',lw='1')
# plt.plot(n/wr1.n0, 135*(wr1.Eneutr(n)/n - C1.M[0]),c='blue')
# plt.plot(n/wr1.n0, 135*(wr2.Eneutr(n)/n - C1.M[0]),c='blue')
# plt.show()


# wr1.reset()
wr1.setDriver()
n, m, r, mg = wr1.stars()
print 'M_NS_MAX =', max(m)
sleep(3)
fNS = wr1.rho[:,0]
contribNS = wr1.getContrib()
nNS = wr1.n
print max(fNS)
# exit()
wr1.C.phi_meson=1
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



# exit()
# wr2.reset(npoints=npoints)


# f1 = wr1.rho[:, 0]
# f2 = wr2.rho[:, 0]
# n = wr1.n/wr1.n0
# print f1-f2
# 
# lines = plt.plot(n, f1, n, f2)
# plt.legend(lines, ['old','new'], loc=0)
# plt.show()
# 
# lines=plt.plot(f1, map(C1.eta_r, f1), f2, map(C2.eta_r, f2))
# plt.legend(lines, ['old','new'], loc=0)
# plt.show()
# 
# wr2 = wr1
# wr2.reset(hyper=1, npoints=npoints, verbose=1, iter=30, timeout=6)
# 
# lines = plt.plot(wr2.n/wr2.n0, wr2.rho[:,0], wr2.n/wr2.n0, wr2.concentrations())
# plt.legend(lines, ['f','p', 's-','s0','s+','x-','x0'], loc=0)
# plt.show()
# 
# fig, ax = plt.subplots(3, 1)
# ax[0].plot(wr2.n/wr2.n0, wr2.rho[:, 0])
# ax[1].plot(wr2.n/wr2.n0, map(wr2.C.eta_r, wr2.rho[:, 0]))
# ax[2].plot(wr2.n/wr2.n0, map(wr2.C.eta_o, wr2.rho[:, 0]))
# plt.show()
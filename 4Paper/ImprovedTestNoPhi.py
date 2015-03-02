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

# TEST OF RHO KIND 6
C2.rho_kind = 6
C2.rho_ld = 0
C1.rho_ld = 0
C2.rho_sat_a = 0 
C2.rho_sat_f1 = 1
C2.rho_sat_f2 = 1
C2.gamma = 0.


C2.beta = 0
C2.rho_sat_val = 0
C2.rho_f = 0.522
rho_f = 0.522
C2.rho_d = 0
C2.rho_width_f = C2.f0
C2.rho_a = 1
C2.rho_width_power = 0
C2.rho_a0 = 0.5
C2.rho_a1 = 2
C2.rho_a2 = 4
C2.rho_a3 = 0
C2.rho_a4 = -0.5

C2.rho_e = 0

C2.rho = 0
C2.drho = 0
C2.d2rho = 0
C2.d2om = 0
C2.dom = 0
C2.gamma2 = 0

C2.rho_tan_a = 0
C2.rho_tan_b = 0
C2.rho_tan_c = 0

C2.beta1 = 0*1.8
C2.beta2 = 0
C2.c1 = 80.

C1.rho = 0
C1.drho = 0
C1.d2rho = 0
C1.d2om = 0
C1.dom = 0
C1.gamma2 = 0

C1.rho_tan_a = 0
C1.rho_tan_b = 0
C1.rho_tan_c = 0

C1.beta1 = 0*1.8
C1.beta2 = 0
C1.c1 = 80.

npoints = 100

frange = np.linspace(0., 1, 1000)

plt.plot(frange, map(C2.eta_r, frange))
plt.show()

lines = []
labels = []
dlist = np.linspace(0, C2.rho_a0, 10)

fig, ax = plt.subplots()

ax.plot(frange, map(C2.eta_r, frange))

# wr2.testPodsiedlowski()

# plt.show()
n_full = 100

wr2.reset(npoints=n_full, iter=200)
vs1 = np.diff(wr2.P)/np.diff(wr2.E)
np1 = wr2.concentrations()[:,1]
fns1 = wr2.rho[:,0]
etar1 = map(C2.eta_r, fns1)
 

C2.rho_ld = 0

f1 = 0.4
 

print C2.rho, C2.drho, C2.d2rho


C2.rho_a_low = 0.04
C2.rho_b_low = -0.03

fsmall = np.linspace(0., 0.27, 20)

def fun_low(x):
    C2.rho_a_low = x[0]
    C2.rho_b_low = x[1]
    res = np.array(map(lambda z: C2.eta_r(z) - C1.eta_r(z), fsmall))
    return np.sum(res**2)

def f(x):
    C2.rho_a = x[0]
    C2.rho_sat_f1 = x[1]
#     C2.rho_tan_a = x[2]
    C2.rho_a0 = x[3]
    C2.rho_a1 = x[4]
    C2.rho_tan_b = x[5]
    f_match = C2.rho_f-0.05
    res =[C2.eta_r(0.27)-1, derivative(C2.eta_r, 0.27, dx=1e-3) - 
            derivative(C1.eta_r, 0.27, dx=1e-3), C2.eta_r(f_match) - C1.eta_r(f_match)]
    res.append(C2.eta_r(0) - C1.eta_r(0))
    res.append(C2.eta_r(0.1) - C1.eta_r(0.1))
#     res.append(derivative(C2.eta_r, C2.rho_f-0.05, dx=1e-3) - derivative(C1.eta_r, C1.rho_f-0.05, dx=1e-3))
    res = np.array(res)
    return np.sum(res*res)


C2.rho_sat_f1 = 0.4
C2.rho_tan_a = 3.
C2.rho_tan_b = 12
C2.rho_a3 = 2
res = optimize.minimize(f, [C2.rho_a, C2.rho_sat_f1, C2.rho_tan_a, C2.rho_a0, C2.rho_a1, C2.rho_tan_b,
                            C2.rho_a2, C2.rho_a4])
print res
print res.x
# exit()
# print optimize.minimize(fun_low, [0,0])

# C2.rho_tan_a = 7.4
# C2.rho_tan_b = -6.6
# C2.rho_tan_c = -3.27


# C2.rho_tan_a = 2
# C2.rho_tan_b = 10
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

fig, ax = plt.subplots(2,2)
ax[0,0].plot(wr2.n[1:]/wr2.n0, vs1, wr2.n[1:]/wr2.n0, vs2)
ax[1,0].plot(wr2.n/wr2.n0, np1, wr2.n/wr2.n0, np2)
ax[0,1].plot(wr2.n/wr2.n0, fns1, wr2.n/wr2.n0, fns2)
ax[1,1].plot(wr2.n/wr2.n0, etar1, wr2.n/wr2.n0, etar2)
print C2.rho_tan_a, C2.rho_tan_b, C2.rho_tan_c
plt.show()

wr2.solve(f0=C2.f0, K0=240., J0=30.)
print wr2.J(), wr2.L()


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
wr1.C.phi_meson=0
wr1.reset(npoints=npoints, timeout=2, hyper=1)
# wr1.setDriver()

fig, ax = plt.subplots(2, 2)
ax[0,0].plot(wr1.n/wr1.n0, wr1.getContrib(), nNS/wr1.n0, contribNS)
ax[1,0].plot(wr1.n/wr1.n0, wr1.concentrations(), wr1.n/wr1.n0, wr1.rho[:,0], nNS/wr1.n0, fNS)
ax[0,1].plot(wr1.rho[:,0], map(wr1.C.eta_r, wr1.rho[:,0]))

plt.show()

n, m, r, mg, mg2 = wr1.stars_crust_hyper()
print max(m)
# exit()

plt.plot(n/wr2.n0, m)
plt.show()


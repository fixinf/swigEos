import eosWrap as eos
from Wrapper import Wrapper
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
# x = linspace(0.0, 3.0, 100)
# l = plot(x, sin(x),x, sin(2*x),x, sin(3*x))
# colors = [c.properties()['color'] for c in l]
# print colors
# show()

C = eos.KVOR_mod()


C.alpha = 0.85
C.omega_a = 12.0
C.rho_a = 0.0
C.phi_f = 0.35
C.phi_a = -3.5
C.omega_f = 0.55
C.d = -10.0
C.phi_gamma = 10.0

f0 = 0.26

wr = Wrapper(C)
for f in linspace(C.f0, f0, 200):
    C.f0 = f
    wr.solve(f0=C.f0)
    

 

C.SetHyperConstants(2)
print C.X_s[0], C.X_s[2], C.X_s[3]
pause(2)
C.phi_meson = 1
print C.f0
print eos.f_eq(array([wr.n0/2, wr.n0/2]), C)
# pause(1)
dx = 1e-4
n_l = np.linspace(0.00001, 5.0, 500)
res = []
rho = [0.0]
feq = []
f = 0

def U_lambda(x, Y, gamma=None):
    f = 0.0
    res = []
    for n in n_l:
        f_eq_arg = array([(1- x)*n*(1-Y), x*n*(1-Y), Y*n])
#         print 'f_eq_arg=',
        f = eos.f_eq(f_eq_arg, C, f)
        mu_arg = np.array([f,(1- x)*n*(1-Y), x*n*(1-Y), Y*n])
        feq.append(f)
        qqq = eos._E(mu_arg, C)/n - C.M[2]
#         qqq = eos.mu(mu_arg, 3, C) - C.M[2]
        res.append(qqq)
    return n_l, np.array(res)*135.0

def U_sigma_m():
    f = 0.0
    res = []
    for n in n_l:
        f_eq_arg = array([0/2, 0/2, 0.0, 0.0, 0.0, 0.0, n/2, n/2])
#         print 'f_eq_arg=',
#         f_eq_arg = array([n/2, n/2, 0.0, 0.0])
        f = eos.f_eq(f_eq_arg, C, f)
        print f
        mu_arg = np.array([f, 0/2, 0/2, 0,0, 0.0, 0.0, 0.0, n/2, n/2])
#         mu_arg = np.array([f, n/2, n/2, 0,0, 0.0])
        feq.append(f)
        qqq = eos.E(mu_arg, C)/n - C.M[7]
#         qqq = eos.mu(mu_arg, 7, C) - C.M[7]
        res.append(qqq)
    print res
    return np.array(res)*135.0

def U_sigma_p():
    f = 0.0
    res = []
    for n in n_l:
        f_eq_arg = array([0.0, 0.0, 0.0, 0.0, 0.0, n])
#         print 'f_eq_arg=',
        f = eos.f_eq(f_eq_arg, C, f)
        mu_arg = np.array([f, 0.0, 0.0, 0,0, 0.0, 0.0, n])
        feq.append(f)
        qqq = eos._E(mu_arg, C)/n - C.M[3]
#         qqq = eos.mu(mu_arg, 4, C) - C.M[3]
        res.append(qqq)
    return np.array(res)*135.0


balb1 = open('balb_lambda_1.dat', 'r')
x_balb = []
y_balb = []
for line in balb1:
    print line
    a,b = line.split()
    x_balb.append(a)
    y_balb.append(b)

x_balb = np.array(x_balb, dtype=float64)*wr.n0/0.16

stoks_L = open('LinL_Stoks', 'r')
x_L_stoks = []
y_L_stoks = []
for line in stoks_L:
    a,b = line.split()
    x_L_stoks.append(a)
    y_L_stoks.append(b)

x_L_stoks = np.array(x_L_stoks, dtype=float64)*wr.n0/0.16

stoks_S = open('XinX_Stoks', 'r')
x_S_stoks = []
y_S_stoks = []
for line in stoks_S:
    a,b = line.split()
    x_S_stoks.append(a)
    y_S_stoks.append(b)

x_S_stoks = np.array(x_S_stoks, dtype=float64)*wr.n0/0.16

print x_balb.dtype
n_l, res = U_lambda(0.0, 0.0)


fig, axx = plt.subplots(1, 2)
ax = axx[0]
plt.subplots_adjust(left=0.15, bottom=0.35)
l, = ax.plot(n_l, res)

def E_frac(_n, x, Y):
    f = 0.0
    res = []
    for n in _n:
        arg_f = np.array([n*(1-Y)*(1-x), n*(1-Y)*x, Y*n])
        f = eos.f_eq(arg_f, wr.C, f)
        arg_e =np.array([f, n*(1-Y)*(1-x), n*(1-Y)*x, Y*n])
        E = eos.E(arg_e, wr.C)/n - C.M[0] 
        res.append(E)
    return np.array(res)*135.0
    
l2, = ax.plot(n_l, E_frac(n_l, 0.5, 0.0))

ax.plot(x_L_stoks, y_L_stoks)
axx[1].plot(x_S_stoks, y_S_stoks)
lSm, = axx[1].plot(n_l, U_sigma_m())
lSp, = axx[1].plot(n_l, U_sigma_p())
axY = plt.axes([0.25, 0.15, 0.35, 0.03])
axX = plt.axes([0.25, 0.25, 0.35, 0.03])
slY = Slider(axY, 'Y', 0.0, 1.0, valinit = 0.2)
slX = Slider(axX, 'x', 0.0, 0.5)
axGamma = plt.axes([0.25, 0.05, 0.35, 0.03])
slGamma = Slider(axGamma, r'$\gamma$', -10.0, 10.0)

axZ = plt.axes([0.65, 0.05, 0.3, 0.03])
slZ = Slider(axZ, 'z', 0.0, 10.0, valinit=0.65)



def Update(val):
    print slX.val, slY.val, slGamma.val
    C.phi_z = slZ.val
    C.phi_gamma = slGamma.val
    n_l, res = U_lambda(slX.val, slY.val, slGamma.val)
#     print len(res)
    l.set_ydata(res)
    res2 = E_frac(n_l, slX.val, slY.val)
    l2.set_ydata(res2)
    lSm.set_ydata(U_sigma_m())
    lSp.set_ydata(U_sigma_p())


slX.on_changed(Update)
slY.on_changed(Update)
slGamma.on_changed(Update)
slZ.on_changed(Update)

plt.show()
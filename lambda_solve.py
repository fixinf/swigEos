import eosWrap as eos
from Wrapper import Wrapper
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy.optimize import optimize
from scipy.optimize import brentq
# x = linspace(0.0, 3.0, 100)
# l = plot(x, sin(x),x, sin(2*x),x, sin(3*x))
# colors = [c.properties()['color'] for c in l]
# print colors
# show()

C = eos.KVOR_mod()


C.alpha = 0.85
C.omega_a = 1.67
C.rho_a = 0.0
C.phi_f = 0.33
C.phi_a = -1.01
C.omega_f = 0.44
C.d = -4.21
C.phi_gamma = 2.0
C.z = 0.65
C.phi_z = 3.5
C.sprime = 0
C.beta = 4.54
C.gamma = 3.78

f0 = 0.26

C.SetHyperConstants(2)

wr = Wrapper(C)
for f in linspace(C.f0, f0, 20):
    C.f0 = f
    wr.solve(f0=C.f0)

print C.eta_o(f0)

C.SetHyperConstants(2)

C.sprime = 1
C.Csp = 400.0
print eos.f_eq(np.array([wr.n0/2.0, wr.n0/2.0]), np.array([0.0, 0.0]), 2, C)
pot =  eos.potentials(np.array([0.26, wr.n0/2, wr.n0/2]), 5, C)
print pot, (C.X_s[2]*pot[0] + C.X_o[2]*pot[2])*135.0
pause(2)

sp = 1 + C.sprime
C.phi_meson = 1



print eos.stepE(0.5, array([0.5]), array([1e-5]), 1, 300, C)
#
dx = 1e-4
n_l = np.linspace(0.00001, 4.0, 1000)
res = []
rho = [0.0]
feq = []

f = np.array([0.0 for j in range(sp)])
flist = []
def U_lambda(n):
    global f 
    res = []
    f_eq_arg = array([0, 0, n])
    f = eos.f_eq(f_eq_arg, f, sp, C)
    flist.append(f)
    mu_arg = np.array([0,0, n])
    mu_arg = np.insert(mu_arg, 0, f)
    feq.append(f)
    qqq = eos._E(mu_arg, C)/n - C.M[2]
#         qqq = eos.mu(mu_arg, 3, C) - C.M[2]
    res.append(qqq)
    return np.array(res)*135.0

print eos.f_eq(np.array([0.0, 0.0, 0.5]), np.array([0.0, 0.0]), 2, C)

l1,=plot(n_l/wr.n0, map(U_lambda, n_l))
plot(n_l/wr.n0, flist)
C.sprime = 0
sp = 1
flist=[]
f = np.array([0.0])
l2,=plot(n_l/wr.n0, map(U_lambda, n_l))
plot(n_l/wr.n0, flist)
# ylim([-6.0, 15.0])
# xlim([0.0, 2.0])

ylabel(r'$U_\Lambda$')
xlabel(r'$n/n_0$')
legend([l1, l2], [r'With $\sigma^*$', r'Without $\sigma^*$'])

show()
C.sprime = 1
sp=2
wr.reset(hyper=1)
print U_lambda(0.5)

# def func_solve(x):
#     C.Csp = x
#     nmin = optimize.fmin(U_lambda, wr.n0, full_output=0)
#     print nmin
#     res = U_lambda(nmin) + 5.0
#     print res
#     return res
# 
# print func_solve(C.Csp)
#     
# x = brentq(func_solve, 100.0, 300.0)
# print x
# 
# print func_solve(x)


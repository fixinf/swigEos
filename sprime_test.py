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
C.SetHyperConstants(2)

C.alpha = 0.85
C.z = 0.65

C.omega_a = 1.67
C.omega_f = 0.44

C.rho_a = 0.0

C.phi_a = -1.01
C.phi_f = 0.33

C.d = -4.21

C.phi_gamma = 3.0
C.phi_z = 3.5

C.sprime = 0
C.Csp = 380.0

C.beta = 4.54
C.gamma = 3.78



f0 = 0.26

nmax = 3.3
C.sprime = 0
wr = Wrapper(C)
for f in linspace(C.f0, f0, 20):
    C.f0 = f
    wr.solve(f0=C.f0)
C.sprime = 1
C.SetHyperConstants(2)

sp = 1 + C.sprime
C.phi_meson = 1

wr.reset(hyper=1, nmax=nmax, npoints=400)
rho = []
sums = []
for r in wr.rho:
    rho.append(r[1+C.sprime:]/sum(r[1+C.sprime:]))
   
for r in rho:
    s = 0.0
    for i, n in enumerate(r):
        s += C.Q[i]*n
    sums.append(s)
flist = wr.rho[:,:1+C.sprime]
mlist = []
masses = []
for f in flist:
    masses = []
    for i in range(8):
        print i
        marg = C.X_s[i]*C.M[0]/C.M[i]*f[0] #+ C.X_sp[i]*C.M[0]/C.M[i]*f[1]
        if C.sprime:
            marg += C.X_sp[i]*C.M[0]/C.M[i]*f[1]
        masses.append(C.phi_n(i,marg))
    mlist.append(array(masses))

mlist = array(mlist)
print mlist.shape
lines = plot(wr.n/wr.n0, mlist)
legend(lines, ['n','p',r'$\Lambda$',r'$\Sigma^-$',r'$\Sigma^0$',
                                         r'$\Sigma^+$', r'$\Xi^-$', r'$\Xi^0$'])
show()

plot(wr.n/wr.n0, rho)
plot(wr.n/wr.n0, wr.rho[:,:1+C.sprime], linestyle='-.')
plot(wr.n/wr.n0, sums, linestyle='--')
show()



print eos.stepE(0.5, array([0.5]), array([1e-5]), 1, 300, C)
#
dx = 1e-4
n_l = np.linspace(0.00001, 3.5, 1000)
res = []
rho = [0.0]
feq = []

init = np.array([0.0 for j in range(sp)])

def U_lambda(x, Y, gamma=None):
    f = init
    res = []
    for n in n_l:
        f_eq_arg = array([(1- x)*n*(1-Y), x*n*(1-Y), Y*n])
#         print 'f_eq_arg=',
        f = eos.f_eq(f_eq_arg, f, sp, C)
        print f
        mu_arg = np.array([(1- x)*n*(1-Y), x*n*(1-Y), Y*n])
        mu_arg = np.insert(mu_arg, 0, f)
        feq.append(f)
        print mu_arg
        qqq = eos._E(mu_arg, C)/n - C.M[2]
#         qqq = eos.mu(mu_arg, 3, C) - C.M[2]
        res.append(qqq)
    return np.array(res)*135.0

def U_XX():
    f = init
    res = []
    for n in n_l:
        f_eq_arg = array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, n, 0.0])
#         print 'f_eq_arg=',
        f = eos.f_eq(f_eq_arg, f, sp, C)
        print f
        mu_arg = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, n, 0.0])
        mu_arg = np.insert(mu_arg, 0, f)
#         feq.append(f)
        print mu_arg
        qqq = eos._E(mu_arg, C)/n - C.M[6]
#         qqq = eos.mu(mu_arg, 3, C) - C.M[2]
        res.append(qqq)
    return np.array(res)*135.0

def U_SmSm():
    f = init
    res = []
    for n in n_l:
        f_eq_arg = array([0.0, 0.0, 0.0, 0.0, 0.0, n])
#         print 'f_eq_arg=',
        f = eos.f_eq(f_eq_arg, f, sp, C)
        print f
        mu_arg = np.array([0.0, 0.0, 0.0, 0.0, 0.0, n])
        mu_arg = np.insert(mu_arg, 0, f)
#         feq.append(f)
        print mu_arg
        qqq = eos._E(mu_arg, C)/n - C.M[4]
#         qqq = eos.mu(mu_arg, 3, C) - C.M[2]
        res.append(qqq)
    return np.array(res)*135.0

def U_LN():
    f = init
    res = []
    for n in n_l:
        f_eq_arg = array([n/2, n/2])
#         print 'f_eq_arg=',
        f = eos.f_eq(f_eq_arg, f, sp, C)
#         print f
        mu_arg = np.array([n/2, n/2, 0.0])
        mu_arg = np.insert(mu_arg, 0, f)
#         feq.append(f)
#         print mu_arg
#         qqq = eos._E(mu_arg, C)/n - C.M[4]
        qqq = eos.mu(mu_arg, sp + 2, C) - C.M[2]
        res.append(qqq)
    return np.array(res)*135.0

# print U_lambda(0.0, 1.0) - U_lambda(0.5, 0.0)
fig, ax = subplots(1, 3)
ax[0].plot(n_l/wr.n0, U_lambda(0.0, 1.0), label='ULL')
ax[0].plot(n_l/wr.n0, U_XX(), label='UXX')
ax[0].plot(n_l/wr.n0, U_SmSm(), label='USMSM')
ax[0].set_ylim([-150, 300.0])
ax[1].plot(n_l/wr.n0, feq)
ax[1].set_ylim([0.0, 2.0])
ax[2].plot(n_l/wr.n0, U_LN())
show()
wr.setDriver()
n_star, m_star, r_star = wr.stars()
print max(m_star)
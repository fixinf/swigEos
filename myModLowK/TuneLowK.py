import matplotlib
from operator import indexOf
matplotlib.use('Qt4Agg')
import numpy as np
import matplotlib.pyplot as plt
from Wrapper import Wrapper
import Models

K0 = 240.
C = Models.myModLowerK(K0)
wr = Wrapper(C)

C.phi_a = -0.
C.phi_f = 0.3

C.omega_f = 0.55
C.omega_a = 0.8
C.omega_c = -20000
C.alpha = 0.4

C.Cs = 234.1580555799
C.Co = 134.8826104616
C.b = 0.0046776700
C.c = -0.0029781609

C.d = -0.5

wr.solve(f0=C.f0, K0=K0)

xu = []
upper = []
with open('klahnUpper', 'r') as f:
    for line in f:
        _x, _u = line.split()
        xu.append(float(_x)/0.16)
        upper.append(float(_u))

xl = []
lower = []
with open('klahnLower', 'r') as f:
    for line in f:
        _x, _l = line.split()
        xl.append(float(_x)/0.16)
        lower.append(float(_l))
        
nrange = np.linspace(0., 4., 400)
P = wr.Psymm(nrange)
E, fs = wr.Esymm(nrange, ret_f=1)
fig, ax = plt.subplots(1, 3)
ax[0].semilogy(nrange/wr.n0, P)
ax[0].semilogy(xu, upper, xl, lower, c='gray')

ax[1].plot(nrange[1:], np.diff(P)/np.diff(E) /wr.const / wr.m_pi**4)
ax[2].plot(nrange, fs)
plt.show()
wr.reset(npoints=1000)
wr.setDriver()
n, m, r, mg = wr.stars()
print max(m), n[np.argmax(m)]/wr.n0
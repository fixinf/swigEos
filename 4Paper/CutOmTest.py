import matplotlib
matplotlib.use('QT4Agg')
import numpy as np
import matplotlib.pyplot as plt
import Models
from Wrapper import Wrapper

C1 = Models.KVOR_cut_mod()
C2 = Models.KVOR_cut_03()
wr1 = Wrapper(C1)
wr2 = Wrapper(C2)

# C1.rho_f = 1.
# C1.gamma = 0.

C2.rho_kind = 2
C2.rho_sat_val = 2.
C2.rho_sat_f1 = 0.53;
C2.rho_sat_f2 = 0.58;
C2.rho_sat_a = 10.

# wr.reset(hyper=1, nmax=2.68)
# fHy = wr.rho[:, 0]
# wr.reset(hyper=0, nmax=2.68)
# fNS = wr.rho[:, 0]
# frange = np.linspace(0., 1, 100)
# 
# fig, ax = plt.subplots(1, 3)
# 
# linesF = ax[0].plot(wr.n/wr.n0, fNS, wr.n/wr.n0, fHy)
# ax[0].legend(linesF, ['NS', 'Hyperons'], loc=0)
# ax[1].plot(frange, map(C.eta_r, frange))
# linesEtaR = ax[2].plot(fNS, map(C.eta_r, fNS), fHy, map(C.eta_r, fHy))
# ax[2].legend(linesEtaR, ['NS', 'Hyperons'], loc=0)
# plt.show()
# 
# fig, ax = plt.subplots(1, 3)
# ax[0].plot(fNS, map(C.eta_o, fNS), fHy, map(C.eta_o, fHy))
# ax[1].plot(fNS, map(C.eta_r, fNS), fHy, map(C.eta_r, fHy))
# ax[2].plot(fNS, map(C.eta_s, fNS), fHy, map(C.eta_s, fHy))
# plt.show()


frange = np.linspace(0., 1., 1000)
plt.plot(frange, map(C1.eta_o, frange))
plt.plot(frange, map(C2.eta_o, frange))
plt.show()

wr1.reset(npoints=100, timeout=2)
wr2.reset(npoints=100)

f1 = wr1.rho[:, 0]
f2 = wr2.rho[:, 0]
n = wr1.n/wr1.n0
print f1-f2

lines = plt.plot(n, f1, n, f2)
plt.legend(lines, ['old','new'], loc=0)
plt.show()

lines=plt.plot(f1, map(C1.eta_o, f1), f2, map(C2.eta_o, f2))
plt.legend(lines, ['old','new'], loc=0)
plt.show()

# wr2 = wr1
wr2.reset(hyper=1, npoints=100, verbose=1, iter=30, timeout=6)

lines = plt.plot(wr2.n/wr2.n0, wr2.rho[:,0], wr2.n/wr2.n0, wr2.concentrations())
plt.legend(lines, ['f','p', 's-','s0','s+','x-','x0'], loc=0)
plt.show()

fig, ax = plt.subplots(3, 1)
ax[0].plot(wr2.n/wr2.n0, wr2.rho[:, 0])
ax[1].plot(wr2.n/wr2.n0, map(wr2.C.eta_r, wr2.rho[:, 0]))
ax[2].plot(wr2.n/wr2.n0, map(wr2.C.eta_o, wr2.rho[:, 0]))
plt.show()
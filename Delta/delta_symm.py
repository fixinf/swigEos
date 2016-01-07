from scipy import optimize
from Wrapper2 import Model

__author__ = 'const'
import Models2
import numpy as np
from matplotlib import pyplot as plt
import Models
wr = Models2.MKVOR_d()
# wr = Models2.Wal_d()
import eosWrap as eos
# C = Models.Cubero
# wr = Model(C)
# plt.plot(wr.sym.nrange, wr.sym.Ebind(wr.sym.nrange))
# plt.show()

m = wr.delta

xs = 1.4
xo = 1.
xr = 1.

mu_dd = []
f=0.
flist = []
elist = []

wr.setDeltaConst(np.array([xs for i in range(4)]),
                 np.array([xo for i in range(4)]),
                 np.array([xr for i in range(4)]),
                 's=%.2f o=%.2f'%(xs, xo))

print([i for i in m.C.X_s])
_n = 2*m.n0
##0.509056100737258
# print(eos.f_eq(np.array([m.n0/2, m.n0/2]), np.array([.2]), 1, m.C)[0])
n_in = [_n/4, _n/4, 0., 0., 0., 0., 0., 0., _n/8, _n/8, _n/8, _n/8]
# _f = eos.f_eq(np.array(n_in), np.array([.6]), 1, m.C)[0]
_f = .6
print(_f)
eparts = np.array([0. for i in range(9)])
print(eos.mu(np.array([_f] + n_in), 10, m.C))
print(eos._E(np.array([.6] + n_in), m.C, eparts))
              # Ef, Eu, Ekn, Ekp, Eo
print(np.sum(eparts) + 4*eos.kineticInt(_n/8, m.C.M[10] - m.C.X_s[10]*m.C.M[0]*_f, 4))
print(eparts) #[ 1.7053506   0.32372391  0.80508056  0.80508056  1.67796202  0.,0., 0.,0.]
# Es, Eo, Eu, Ekn, Ekd
# 1.7053505968483933 1.677962022110554 0.3237239070247516 1.6101611206767723 1.743659225476203
#Различаются только E_kin(Delta) ???
print(eos.kineticInt(_n/4, m.C.M[0] - m.C.X_s[0]*m.C.M[0]*_f, 2))
print(4*eos.kineticInt(_n/8, m.C.M[10] - m.C.X_s[10]*m.C.M[0]*_f, 4)) #1.743659225476203



# And where's the difference? o0
# exit()
# for n in m.nrange:
#     f = eos.f_eq(np.array([n/4, n/4, 0., 0., 0., 0., 0., 0., n/8, n/8, n/8, n/8]), np.array([f]), 1, m.C)[0]
#     flist.append(f)
#     mu_dd.append(m.m_pi*eos.mu(np.array([f, n/4, n/4, 0., 0., 0., 0., 0., 0., n/8, n/8, n/8, n/8]), 10, m.C))
#     elist.append(eos._E(np.array([f, n/4, n/4, 0., 0., 0., 0., 0., 0., n/8, n/8, n/8, n/8]), m.C))
#
# np.savetxt("mu_test.dat", np.array([m.nrange/m.n0, mu_dd, flist, elist]).transpose(), fmt='%.6f')
# exit()
# frac, E, Ebind, Eparts, fd = m.dumpDeltaSym()
# nd = frac[:, -1]
#
# plt.plot(m.nrange, fd)
# plt.plot(m.nrange, Eparts[:,1])
# plt.show()




# for xs, xo in [[1.25, 1.], [1., 1.], [1.15, 0.9]]:
for xs, xo in [[1.4, 1.]]:
    wr.setDeltaConst(np.array([xs for i in range(4)]),
                     np.array([xo for i in range(4)]),
                     np.array([1., 1., 1., 1.]),
                     's=%.2f o=%.2f'%(xs, xo))
    m.dumpDeltaSym('Ebind')
    m.inspect_f()
from scipy.misc.common import derivative

__author__ = 'const'
from Nonlin.nonlin import NLVector
import numpy as N

# NL3
# gs = 104.3871
# gv = 165.5854
# gr = 79.6
# kappa = 3.8599
# lamb = -0.015905
# mpi = 135
# mn = 939/ mpi
# mo = 782.501/ mpi
# mr = 763/ mpi
# ms = 508.194 / mpi
# m = NLVector(mn**2 * gs / ms**2, mn**2 * gv / mo**2, mn**2 * gr / mr**2, kappa/2/mn/mpi, lamb/6, 0, 0)

# FSU

mpi = 135
ms = 491.5 / mpi
mo = 782.5 / mpi
mr = 763 / mpi
gs = 112.1996
gv = 204.5469
gr = 138.4701
kappa = 1.4203
lamb = 0.023762
zeta = 0.06
L_v = 0.030
mn = 939 / mpi

m = NLVector(mn**2 * gs / ms**2, mn**2 * gv / mo**2, mn**2 * gr / mr**2,
             kappa/2/mn/mpi, lamb/6, zeta*gv**2/24, L_v*gr*gv)

m2 = NLVector(mn**2 * gs / ms**2, mn**2 * gv / mo**2, mn**2 * gr / mr**2,
             kappa/2/mn/mpi, lamb/6, zeta*gv**2/24, L_v*gr*gv, cut_a=1., cut_c=0.3)

m = NLVector(196.3428132661, 90.7681870142, 88.7261140316, 0.0089455122, 0.0077076578, 0, 0)
m.n0 = 0.16*(197.33/135)**3
m.mn = 938./135
#
# gs = 10.217
# go = 12.868
# gr = 4.474
# mpi = 135
# mn = 939/ mpi
# mo = 782.501/ mpi
# mr = 763/ mpi
# ms = 508.194 / mpi
# g2 = -10.431 * (197.33 / 939 / gs**3)
# g3 = -28.885 / gs**4
#
# print(mn**2 * gs**2 / ms**2, mn**2 * go**2 / mo**2, mn**2 * gr**2 / mr**2, g2, g3)
# m = NLVector(mn**2 * gs**2 / ms**2, mn**2 * go**2 / mo**2, mn**2 * gr**2 / mr**2, g2, g3, 0, 0)
# m = NLVector(373.176, 245.458, 149.67, -0.24578e-2, -0.34334e-2, 0, 0)



n0 = m.n0
print(m.f_eq(m.n0/2, m.n0/2))
print(m.vec_eq(m.n0/2, m.n0/2))
# print(m.go * n0 / m.mo**2)
nrange = N.linspace(0, 4, 100)
print(135*(m.E(m.n0/2, m.n0/2)/m.n0 - m.mn))
print(derivative(lambda z: m.E(z/2, z/2)/z, n0, dx=1e-3, order=7))
print(m.K(m.n0))
# print(m.J(2*m.n0))
print(m.J(m.n0))
# n_p, n_e, n_mu = m.np_eq(0.47/0.148 * n0)
# DU = 1 / (1 + (1 + (n_e/(n_e + n_mu))**(1/3))**3)
nnp = 2
print(N.array(m.np_eq(nnp*n0)) / (nnp*n0))
# print(m2.K(m.n0))
# print(m2.J(m.n0))
# plt.plot(nrange/m.n0, list(map(lambda z: 135*(m.E(z/2, z/2)/z - m.mn), nrange)))
# plt.plot(nrange/m.n0, list(map(lambda z: 135*(m2.E(z/2, z/2)/z - m.mn), nrange)))
# plt.show()
# print(m.P(nrange, N.array(list(map(lambda z: m.E(z/2, z/2), nrange)))))




__author__ = 'const'
from Nonlin.nonlin import NLVector
import numpy as N
from matplotlib import pyplot as plt

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
#
# m = NLVector(196.3428132661, 90.7681870142, 88.7261140316, 0.0089455122, 0.0077076578, 0, 0)
# m.tabMassCrust('FSUgold')


# exit()


m2 = NLVector(mn**2 * gs / ms**2, mn**2 * gv / mo**2, mn**2 * gr / mr**2,
             kappa/2/mn/mpi, lamb/6, zeta*gv**2/24, L_v*gr*gv, cut_a=1., cut_c=0.4)

m2.tabEos(name='FSUgold_cut04', npoints=100)

exit()

mpi_2_mevfm3 = 135 * (135 / 197.33)**3

clist = [1.0, 0.2, 0.3, 0.4]


E, P, n, np = m.EPN_NS(50, ret_np=1)
E *= mpi_2_mevfm3
P *= mpi_2_mevfm3
n /= m.n0

N.savetxt('test.dat', N.array([n, E, P, np]).transpose())



plt.plot(n/m.n0, np)
plt.show()
plt.plot(E, P)
plt.show()

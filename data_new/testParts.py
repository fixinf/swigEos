__author__ = 'const'
import matplotlib
matplotlib.use('QT4Agg')
import Models2
import time

from matplotlib import pyplot as plt
import numpy as np
m = Models2.KVOR()
wr = m.sym
t1 = time.time()
E, P, n = wr.EPN()
C = wr.C
t2 = time.time()
print(t2 - t1)
# exit()
Elist, flist = wr.E(nrange=m.nrange, ret_f=1)
Eparts = wr.Eparts
Pparts = wr.Pparts

part_lines = plt.plot(n, Pparts, linestyle='--', lw='3')
plt.plot(n, np.sum(Pparts, axis=1).transpose(), label='sum')
plt.plot(n, P, label='True')
plt.legend()
plt.legend(part_lines, ['f','U','K_n', 'K_p', 'om', 'rho'])

plt.show()
print(P/wr.mpi4_2_mevfm3 - np.sum(Pparts, axis=1).transpose())
# plt.plot(n, flist)
# plt.show()
plt.plot(n, -C.M[0]**4/ 2 / C.Cs * flist**2 * np.array(list(map(C.eta_s, flist))), label='ppanal')
plt.plot(n, C.Co * n**2 / (2 * C.M[0]**2) / np.array(list(map(C.eta_o, flist))), label='om_anal')
# print(Pparts[:,0])
plt.plot(n, Pparts[:, 0], label='ppnum')
plt.plot(n, Pparts[:, 4], label='om_num')
plt.legend()
plt.show()
exit()
elines = plt.plot(n, Eparts, linestyle='--', lw=3)
plt.legend(elines, ['f','U','K_n', 'K_p', 'om', 'rho'])
# plt.plot(n, E)



print(Eparts)
print(E - np.sum(Eparts, axis=1))
plt.plot(n, np.sum(Eparts, axis=1).transpose())
plt.show()
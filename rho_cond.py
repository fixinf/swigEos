__author__ = 'const'
import matplotlib
matplotlib.use("QT5Agg")
from rho_wrap import *
import numpy as np
import eosWrap as eos


from matplotlib import pyplot as plt
import Models2

wr = Models2.myMod()
C = wr.nucl.C

nrange = np.linspace(0.5, 2.5, 1000)
init = [0.01, 0.2, 0.5]
nr_list = []
np_list = []
nch_list = []
flist = []
for n in nrange:
    res = eq_sol(n, C, init=init)
    init = res[0]
    for i, n in enumerate(init):
        if n < 0:
            init[i] = 1e-6
#     print(res)
    np_list.append(init[0])
    nr_list.append(res[1][0])
    nch_list.append(res[1][1])
    flist.append(init[1])
print(nr_list)
print(nch_list)

plt.plot(nrange/C.n0, np.array(nr_list)/ C.n0, label='nr_list')
plt.plot(nrange/C.n0, (nrange - 2 * np.array(np_list))/C.n0, label='nn-np')
plt.plot(nrange/C.n0, np.array(nch_list)/C.n0, label='n_ch')
plt.plot(nrange/C.n0, np.array(flist), label='f')

# wr.nucl.reset()
# plt.plot(wr.nrange/wr.n0, wr.nucl.rho[:,0])
# plt.plot(wr.nrange/wr.n0, (wr.nucl.rho[:,1] - wr.nucl.rho[:,2])/wr.n0)
# plt.ylim([0., 10.])
plt.legend()
plt.show()
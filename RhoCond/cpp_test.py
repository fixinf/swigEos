import eosWrap as eos
import numpy as np
from matplotlib import pyplot as plt

import Models2

wr = Models2.KVOR()
C = wr.C

init = np.array([0.22042071])
finit = np.array([0.1])

init2 = np.array([0.22042071, 1.57862672])

nrange = np.linspace(0*wr.n0, 6*wr.n0, 1000)
# nrange = [1.8]
arr = []
res = init2
finit = np.array([0.43540519])

for n in nrange:
    print(n/wr.n0)
    res = eos.stepE_rho(n, init2, finit, len(init2)+1, 100, 0.1, C)
    res2 = eos.stepE(n, init, finit, len(init), 100, C)
    init2 = np.array(res[:2])
    nn = n - res[0]
    n_in_f = np.array([nn, res[0]])
    finit = eos.f_eq_rho(n_in_f, np.array(finit), 1, res[1], C)
    arr.append(res)
    # res2 = eos.stepE(n, init, finit, len(init), 30, C)
    print(res2)
    print(res, finit)

    # print(res2)

arr = np.array(arr)
plt.plot(nrange/wr.n0, arr[:, 2] / wr.n0)
plt.show()

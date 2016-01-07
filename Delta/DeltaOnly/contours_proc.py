import eosWrap as eos
import Models2
import numpy as np
from matplotlib import pyplot as plt

import joblib as jl

wr = Models2.Wal_d()
m = wr.delta_phi

def fun(xs, xo):
    wr.setDeltaConst(np.array([xs for i in range(4)]),
                 np.array([xo for i in range(4)]),
                 np.array([1., 1., 1., 1.]),
                 's = %.2f o = %.2f' % (xs, xo))
    print(m.foldername)
    # exit()
    m.loadEos()
    lines = plt.plot(m.nrange/m.n0, m.concentrations())
    # print(m.concentrations())
    # exit()
    # plt.legend(lines, m.part_names)
    # plt.show()
    ##### Pressure behavior
    trans = 0
    conc = m.concentrations()
    E,P,N = m.EPN()
    for i, p in enumerate(P):
        if i > 0:
            if P[i-1] > P[i]:
                trans = 1

    ncrit = np.zeros(conc.shape[1])
    ncrit[:] = 8.
    for n, col in enumerate(conc.transpose()):
        for i, frac in enumerate(col):
            if frac > 1e-6:
                ncrit[n] = m.nrange[i]/m.n0
                break
    return [trans, ncrit]

npoints = 10
out = []
xslist = np.linspace(0.8, 1.3, npoints)
xolist = np.linspace(0.8, 1.3, npoints)

for xs in xslist:
    for xo in xolist:
        print(xs, xo)
        res = fun(xs, xo)
        out.append(np.insert(res[1], 0, [xs, xo, res[0]]))

np.savetxt("test_100.dat", np.array(out))
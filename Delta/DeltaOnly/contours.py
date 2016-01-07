__author__ = 'const'
import Models2
from matplotlib import pyplot as plt
import numpy as np
import eosWrap as eos
import joblib as jl
wr = Models2.Wal_d()
m = wr.delta_only
wr.dumpProps()
# exit()
def fun(xs, xo):
    wr.setDeltaConst(np.array([xs for i in range(4)]),
                 np.array([xo for i in range(4)]),
                 np.array([1., 1., 1., 1.]),
                 's = %.2f o = %.2f' % (xs, xo))
    print(m.foldername)

    # exit()
    m.reset(timeout=None, iterations=60)
    m.dumpEos()
    m.dumpMu()
    m.dumpMeff()
    # m.dumpScalings()
    m.dumpVs()
    m.dumpParts()
    ##### Pressure behavior
    trans = 0
    conc = m.concentrations()
    E,P,N = m.EPN()
    for i, p in enumerate(P):
        if i > 0:
            if P[i-1] > P[i]:
                trans = 1

    ncrit = np.zeros(conc.shape[1])
    ncrit[:] = 10.
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
res = jl.Parallel(n_jobs=4)(jl.delayed(fun)(i,j) for i in xslist for j in xolist)
exit()
for xs in xslist:
    for xo in xolist:
        print(xs, xo)
        res = fun(xs, xo)
        out.append(np.insert(res[1], 0, [xs, xo, res[0]]))

np.savetxt("test_4.dat", np.array(out))
__author__ = 'const'
from Wrapper2 import Model

__author__ = 'const'
import Models2
from matplotlib import pyplot as plt
import numpy as np
import eosWrap as eos

wr1 = Models2.myMod()
wr2 = Models2.MKVORexp()
# wr1 = Model(eos.KVOR_d)
m1 = wr1.hyper
m2 = wr1.delta_phi

wr = Models2.MKVOR_d()
m = wr.delta_phi

def fun(xs, xo):
    wr.setDeltaConst(np.array([xs for i in range(4)]),
                 np.array([xo for i in range(4)]),
                 np.array([1., 1., 1., 1.]),
                 'o')
    m.reset()
    ##### Pressure behavior
    trans = 0
    conc = m.concentrations()
    E,P,N = m.EPN()
    for i, p in enumerate(P):
        if i > 0:
            if P[i-1] > P[i]:
                trans = 1

    ncrit = np.zeros(conc.shape[1])
    for n, col in enumerate(conc.transpose()):
        for i, frac in enumerate(col):
            if frac > 1e-6:
                ncrit[n] = m.nrange[i]/m.n0
                break
    return [trans, ncrit]

# print(fun(1.147368421052631593e+00, 8.947368421052630527e-01))
#
# exit()

npoints = 10
out = []
for xs in np.linspace(0.8, 1.2, npoints):
    for xo in np.linspace(8.947368421052630527e-01, 8.947368421052630527e-01, npoints):
        print(xs, xo)
        res = fun(xs, xo)
        out.append(np.insert(res[1], 0, [xs, xo, res[0]]))

np.savetxt("test_1.dat", np.array(out))



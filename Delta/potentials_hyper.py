import eosWrap as eos
from os.path import join
from scipy import optimize
from scipy.interpolate.interpolate import interp1d
import Models2
import numpy as np
from matplotlib import pyplot as plt

import joblib as jl

c = 0.3
# wr = Models2.Cubero_cut(c, suffix="c=%.1f _Delta"%c)
# wr = Models2.KVORcut03(J0=30.)
# wr = Models2.Wal_d(J0=30., f0=0.3)
# wr = Models2.KVOR_d()
wr = Models2.MKVOR_d()
m = wr.delta
C = m.C
S,V = wr.dumpPotentials()
wr.dumpProps()
DUs = []
def proc(U):
    i = 8
    iS = interp1d(m.nrange/m.n0, S)
    iV = interp1d(m.nrange/m.n0, V)
    xs_d = (U - C.X_o[i]*iV(1.))/iS(1.)
    print(xs_d)

    xo=1.

    wr.setDeltaConst(np.array([xs_d for i in range(4)]),
                 np.array([xo for i in range(4)]),
                 np.array([1., 1., 1., 1.]),
                 's = %.2f o = %.2f U = %.0f' % (xs_d, xo, U))


    wr.delta_phi.loadEos()
    wr.delta_phi.dumpMu()
    wr.delta_phi.dumpMeff()
    wr.delta_phi.dumpMassesCrust()

    # wr.dumpPf()
    wr.delta_phi.dumpScalings()

ulist = [-100., -150.]
# ulist = [0.]
# res = jl.Parallel(n_jobs=2)(jl.delayed(proc)(u) for u in ulist)
res = [proc(u) for u in ulist]
# with open(join(wr.foldername, 'DUs.dat'), 'w') as f:
#     for i in DUs:
#         f.write(str(i)+'\n')
#
#
#

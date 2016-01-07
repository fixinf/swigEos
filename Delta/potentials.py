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
wr = Models2.KVORcut03()
# wr = Models2.Wal_d(J0=28.)
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
    print(iV(1.), iS(1.))

    xs_d = (U - C.X_o[i]*iV(1.))/iS(1.)
    print(xs_d)
    # exit()
    xo=1.

    wr.setDeltaConst(np.array([xs_d for i in range(4)]),
                 np.array([xo for i in range(4)]),
                 np.array([1., 1., 1., 1.]),
                 's = %.2f o = %.2f U = %.0f' % (xs_d, xo, U))

    # wr.delta_sym.dumpEos()
    # wr.delta_only.dumpEos()
    # wr.delta_only.dumpMu()
    # wr.delta_only.dumpMeff()
    # wr.delta_sym.dumpEos()
    # wr.delta_only.dumpEos()
    # wr.delta_only.dumpMu()
    # wr.delta_only.dumpMeff()
    wr.delta_phi.loadEos()
    wr.delta_phi.dumpMassesCrust()
    exit()
    pfs = wr.delta_only.dumpPf(write=1)
    i_start = 10
    dpf = pfs[i_start:, 0] - pfs[i_start:, 1] - pfs[i_start:, -2]
    # lines = plt.plot(wr.nrange[i_start:]/wr.n0, np.abs(dpf))
    # plt.show()
    # print(np.argmin(np.abs(dpf)))
    # print(m.nrange[350])
    # print(m.nrange[350]/m.n0)
    index = 0
    for i, elem in enumerate(dpf):
        if elem < 0:
            index = i
            break

    n_crit = m.nrange[index]

    iDU = interp1d(wr.nrange[i_start:], dpf)
    n_crit2 = optimize.bisect(iDU, n_crit-0.5, n_crit+0.5)


    Mdu = 0.
    try:
        masses = np.loadtxt(join(wr.delta_only.foldername, wr.delta_only.filenames['mass_crust']), skiprows=1)
        iM = interp1d(masses[:, 0], masses[:, 1])
        Mdu = iM(n_crit2/wr.n0)
    except FileNotFoundError:
        pass



    print(wr.delta_only.foldername, n_crit/m.n0, n_crit2/m.n0, Mdu)
    DUs.append([wr.delta_only.foldername, n_crit/m.n0, n_crit2/m.n0, Mdu])
    # print(wr.delta_only.xDUp)
    # exit()
    # wr.delta_sym.dumpEos()
    # print(wr.delta_phi.foldername)
    # try:
    #     # wr.delta_phi.loadEos()
    # except FileNotFoundError:
    # wr.delta_phi.dumpEos()
    # ne, nmu, dus = wr.delta_only.lepton_concentrations(ret_du=1)
    # np.savetxt(join(wr.delta_only.foldername, 'du.dat'), np.array([wr.nrange/wr.n0, dus]).transpose(),
    #            fmt='%.6f')
    # wr.delta_only.dumpMassesCrust()

# ulist = [-110., -120., -135., -140., 0., -90.,  -130., -150.]
ulist = [-50]
# ulist = [0.]
res = jl.Parallel(n_jobs=2)(jl.delayed(proc)(u) for u in ulist)
exit()
res = [proc(u) for u in ulist]
with open(join(wr.foldername, 'DUs.dat'), 'w') as f:
    for i in DUs:
        f.write(str(i)+'\n')




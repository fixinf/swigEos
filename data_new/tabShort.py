__author__ = 'const'
import Models2
from matplotlib import pyplot as plt
import numpy as np
# m = Models2.myMod_L(1.2, 7.5)
models = [Models2.KVORcut02(), Models2.KVORcut04()]
# models = [Models2.KVOR()]#, Models2.myMod()]
for m in models:
    m.dumpEos()
    for n in m.children[:-1]:
        n.dumpEos()
    # m2 = Models2.KVOR()

    # m.dumpUofE()
    # m.dumpUofE_anti()
    # m.sym.dumpJ()
    # # # m.nucl.dumpMassesCrust(nmin=2*m.n0, nmax=3*m.n0, fname='mass_crust_detail.dat')
    # # # # exit()
    # m.sym.dumpScalings()
    # m.nucl.dumpScalings()
    # frange = np.linspace(0., 1., 100)
    # plt.plot(frange, map(m.C.eta_o, frange))
    # plt.plot(frange, map(m2.C.eta_o, frange))
    # plt.show()
    # m.dumpProps()
    # m.dumpScalings()
    # m.dumpEos()
    # m.nucl.dumpMassesCrust()

    # m.dumpAll()
    # m.sym.dumpVs()
    # m.nucl.dumpVs()

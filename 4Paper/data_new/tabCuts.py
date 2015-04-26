__author__ = 'const'
import Models2
from matplotlib import pyplot as plt
import numpy as np
# m = Models2.myMod_L(1.2, 7.5)
K0 = 200.
f0 = 0.12
J0 = 30.
zlist = [1., 0.5, 0.4, 0.3, 0.2]

models = [Models2.Cubero_cut(z, K0=K0, f0=f0, J0=J0, suffix='K=%3.0f=%1.2f'%(K0,f0))
          for z in zlist]
frange = np.linspace(0, 1, 100)
lines = []
legend = []
for m in models:
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
    print m.foldername
    m.dumpProps()
    m.dumpScalings()
    # l, = plt.plot(frange, map(m.C.U, frange))
    # lines.append(l)
    # legend.append('%1.1f'%m.C.c_sigma)
    # m.dumpEos()
    # m.nucl.dumpMassesCrust()
    m.dumpAll(hyper=0)
    # m.sym.dumpVs()
    # m.nucl.dumpVs()
    # m.neutr.dumpVs()
# plt.legend(lines, legend)
# plt.show()
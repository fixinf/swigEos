from scipy import optimize
from scipy.optimize._root import root
from scipy.optimize.zeros import bisect

__author__ = 'const'
import Models2
from matplotlib import pyplot as plt
import numpy as np
plist = [[0.8, 7.5],
         [1.2, 7.5],
         [1.5, 7.5],
         [1.4, 9.]
         ]

def fun_L(p):
    m = Models2.myMod_L(*p)
    print(m.sym.L())
    return m.sym.L()

plist = [[0.,0.], [0.8, 3.6], [0.8, 6]]
# p = [[0., 0.]]
# pp = []
# for _p in p:
#     pp.append([_p, fun_L(_p)])
# print(pp)
# exit()

for p in plist:
    m = Models2.myMod_L(*p)
    m.nucl.dumpMassesCrust(nmin=.6)

    m.dumpEos()
    # E, P, N = m.nucl.EPN()
    # plt.plot(E, P)
    # plt.show()
    # m.nucl.dumpProfiles(2.7*m.n0)
    # m.nucl.dumpMassesCrust(nmin=.6)
    m.dumpUofE()
    m.dumpUofE_anti()
    m.sym.dumpJ()
    # m.nucl.dumpMassesCrust(nmin=2*m.n0, nmax=3*m.n0, fname='mass_crust_detail.dat')
    # # # # # # exit()
    # m.dumpScalings()
    m.sym.dumpScalings()
    m.nucl.dumpScalings()
    # # # frange = np.linspace(0., 1., 100)
    # # # plt.plot(frange, map(m.C.eta_o, frange))
    # # # plt.plot(frange, map(m2.C.eta_o, frange))
    # # # plt.show()
    m.dumpProps()
    # m.dumpEos()
    m.dumpEos()
    # m.dumpAll()
    m.sym.dumpVs()
    m.nucl.dumpVs()
    m.neutr.dumpVs()
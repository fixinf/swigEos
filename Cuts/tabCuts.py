__author__ = 'const'
import Models2
from matplotlib import pyplot as plt
import numpy as np
# m = Models2.myMod_L(1.2, 7.5)
K0 = 250.
f0 = 0.2
J0 = 30.
zlist = [1.]#, 0.5, 0.4, 0.3, 0.2]

params = [
    {
        'K' : 250.,
        'f0' : 0.2,
        'Cs' : 196.3428132661,
        'Co' : 90.7681870142,
        'Cr' : 88.7261140316,
        'b' : 0.0089455122,
        'c' : 0.0077076578
    },

# K :200., f0=0.2
    {
        'K' : 200.,
        'f0' : 0.2,
        'Cs' : 213.2607824762,
        'Co' : 90.7681524182,
        'Cr' : 88.7261190772,
        'b' : 0.0150061019,
        'c' : -0.0124942770
    },
#K=200., f0=0.19
    # {
    #     'K' : 200.,
    #     'f0' : 0.19,
    #     'Cs' : 205.0913029877,
    #     'Co' : 84.4290026333,
    #     'Cr' : 89.6408647791,
    #     'b' : 0.0166215109,
    #     'c' : -0.0114022252
    # },
# K=200., f0=0.3
    {
        'K' : 200.,
        'f0' : 0.3,
        'Cs' : 273.0836626079,
        'Co' : 153.6303271470,
        'Cr' : 78.4149744148,
        'b' : 0.0057865796,
        'c' : -0.0072525183
},
#K=200., f0=0.13
    {
        'K' : 200.,
        'f0' : 0.13,
        'Cs' : 143.7534593933,
        'Co' : 39.8340610093,
        'Cr' : 95.5473126488,
        'b' : 0.0102107661,
        'c' : 0.2689864144
    }
]
for p in params:
    models = [Models2.Cubero_cut(z, K0=p['K'], f0=p['f0'], J0=J0, suffix='K=%3.0ff=%1.2f'%(p['K'], p['f0']), params=p)
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
        print(m.foldername)
        # m.dumpProps()
        # m.dumpScalings()
        # l, = plt.plot(frange, map(m.C.U, frange))
        # lines.append(l)
        # legend.append('%1.1f'%m.C.c_sigma)
        # m.dumpEos()
        # m.nucl.dumpMassesCrust()
        # m.dumpAll(hyper=0)
        # m.nucl.dumpMasses()
        # m.sym.dumpVs()
        # m.sym.dumpParts()
        # m.neutr.dumpParts()
        # m.nucl.dumpParts()
        # m.nucl.dumpVs()
        # m.neutr.dumpVs()
    # plt.legend(lines, legend)
    # plt.show()
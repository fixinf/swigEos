__author__ = 'const'
import Models2
import numpy as np
from matplotlib import pyplot as plt
import eosWrap as eos

models = [Models2.KVOR(), Models2.myMod(), Models2.KVORcut02(),
          Models2.KVORcut03(), Models2.KVORcut04()]

# models = [Models2.KVOR()]
models = [Models2.myMod()]
# models = [Models2.Cubero_cut(0.3)]
J0 = 30.
p =     {
        'K' : 250.,
        'f0' : 0.2,
        'Cs' : 196.3428132661,
        'Co' : 90.7681870142,
        'Cr' : 88.7261140316,
        'b' : 0.0089455122,
        'c' : 0.0077076578
    }



# models=[Models2.Cubero_cut(.2, K0=p['K'], f0=p['f0'], J0=J0, suffix='ololo', params=p)]
# models = [Models2.myMod()]
for m in models:
    # m.hyper_phi_sigma.dumpMu()
    m.nucl.stars_crust(ncut_crust=.45, ncut_eos=.6, npoints=1)
    m.nucl.stars_crust(ncut_crust=.45, ncut_eos=.6, npoints=1)
    # n_p = np.linspace(0, m.n0, 20)
    # print(m.sym.J(), m.sym.K(), m.sym.L(), m.sym.Kprime(), m.sym.Ksymm())
    # print(m.sym.P(n_p)[-1])
    # print(m.sym.Ebind(n_p)[-1])
    # # m.setParams(179.56, 87.600, 100.64, 7.7346e-3, 3.4462e-4, 0.195)
    # m.setParams(234.15, 134.88, 81.842, 4.6750e-3, -2.9742e-3, 0.27)
    # print(m.sym.J(), m.sym.K(), m.sym.L(), m.sym.Kprime(), m.sym.Ksymm())
    # print(m.sym.P(n_p)[-1])
    # print(m.sym.Ebind(n_p)[-1])
    # m.dumpPotentials()
    # m.dumpAll(hyper=0)
    # m.nucl.dumpMassesCrust()
    # m.hyper_phi.dumpEos()
    # m.hyper_phi_sigma.dumpEos()
    # nrange = np.linspace(15.9*m.n0, 16.1*m.n0, 2000)
    # E, f = m.sym.E(nrange, ret_f=1, f=.36)
    # fig, ax = plt.subplots(2, 1)
    # ax[0].plot(nrange/m.n0, f)
    # ax[1].plot(nrange/m.n0, np.gradient(f, nrange[1]-nrange[0]))
    # ax[1].set_ylim([0, 0.2])
    # plt.show()
    # for n_c in [.001, .1, .2, .3, .4, .6, .8, 1., 1.2, 1.4, 1.6, 1.8]:
    # # n, mass, r, mg1, mg2 = m.nucl.dumpMassesCrust(write=0, ret_frac=0, npoints=100)
    # # m.nucl.dumpProfiles(n[np.argmax(mass)]*m.nucl.n0)
    # # m.dumpAll()
    #     n, mg, r, mb1, mb2 = m.nucl.stars_crust(ncut_crust=n_c, ncut_eos=n_c+0.2, inter='linear', nmin=.4)
    #     plt.plot(r, mg, label='%.1f'%n_c)
    #
    # plt.legend()
    # # plt.savefig('masses_crust.jpg')
    # plt.show()
    # for child in m.children:
    #     child.dumpVs()



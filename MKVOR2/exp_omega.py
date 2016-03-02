from cmath import cosh

from scipy.interpolate.interpolate import interp1d
from scipy.optimize.minpack import leastsq

import Models2
import Models
from matplotlib import pyplot as plt
import numpy as np


for wr in [Models2.MKVOR_d(), Models2.MKValpha00(0.68), Models2.myModExpOmega(0.68)]:
    for f in [0.68]:#np.linspace(0.5, 1., 10., endpoint=0):
        # wr = Models2.myModExpOmega(f)
        # wr = Models2.OmegaWide(f)
        # wr.dumpAll(hyper=0)
        # wr.hyper.dumpEos()
        # wr.hyper.dumpEos()
        # exit()
        m = wr.delta_sym
        S,V = wr.dumpPotentials()
        for U in [-30., -50.,  -100.]:
            i = 8
            iS = interp1d(m.nrange/m.n0, S)
            iV = interp1d(m.nrange/m.n0, V)
            print(iV(1.), iS(1.))

            xs_d = (U - m.C.X_o[i]*iV(1.))/iS(1.)
            print(xs_d)
            # exit()
            xo=1.

            wr.setDeltaConst(np.array([xs_d for i in range(4)]),
                         np.array([xo for i in range(4)]),
                         np.array([1., 1., 1., 1.]),
                         's = %.2f o = %.2f U = %.0f' % (xs_d, xo, U))

            m.dumpEos()
            wr.delta_phi.dumpEos()
            wr.delta_phi.dumpMassesCrust()

            wr.delta_phi_sigma.dumpEos()
            wr.delta_phi_sigma.dumpMassesCrust()


            # wr.hyper.dumpEos()
            # frange = np.linspace(0., 1, 100)
            #
            # plt.semilogy(wr.nrange/wr.n0, wr.sym.P(wr.nrange))
            # plt.semilogy(wr.nrange/wr.n0, wr2.sym.P(wr.nrange))
            # plt.ylim([1., 2000.])
            # plt.xlim([1., 8.])
            # plt.show()
            # wr2.dumpEos()
            # wr2.hyper.dumpEos()
            #
            # m.reset()
            # plt.plot(m.nrange[:m.rho.shape[0]]/m.n0, m.concentrations())
            # plt.show()
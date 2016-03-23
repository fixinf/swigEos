import eosWrap as eos
import numpy as np
from matplotlib import pyplot as plt

import Models2

for f in np.linspace(0.6, 0.66, 10):
    wr = Models2.MKValpha00(f)
    m = wr.rcond_hyper_phi
    m.reset()
    m.dumpEos()
    # lc = m.lepton_concentrations().transpose()
    # nc = m.nc / m.nrange
    #
    # lines = plt.plot(m.nrange/m.n0, m.concentrations())
    # plt.legend(lines, m.part_names)
    # plt.plot(m.nrange/m.n0, nc)
    # plt.plot(m.nrange/m.n0, lc)
    #
    # print((m.lepton_concentrations()[0] * m.nrange)[-1])
    # print(m.nc[-1])
    #
    # print(m.concentrations()[:, 1] + m.concentrations()[:, 5] - lc[:, 0] - lc[:, 1] - nc)
    #
    #
    # plt.show()
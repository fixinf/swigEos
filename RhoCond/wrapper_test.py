import eosWrap as eos
import numpy as np
from matplotlib import pyplot as plt

import Models2

wr = Models2.MKValpha00(0.66)
wr.nucl.dumpMassesCrust()
exit()
# wr.dumpAll()
# wr.delta_sym.dumpEos()
# exit()
# wr = Models2.KVOR()
m = wr.rcond_hyper_phi
m.reset()
m.dumpEos()
# m.dumpMassesCrust()
line_nc = plt.plot(m.nrange/wr.n0, m.nc/m.nrange)
lines = plt.plot(m.nrange/m.n0, m.concentrations())
line_f = plt.plot(m.nrange/m.n0, m.rho[:, 0])
plt.legend(line_f + line_nc + lines, ['f'] + ['nc'] + m.part_names)
plt.show()
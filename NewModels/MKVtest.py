import eosWrap as eos
import numpy as np
from matplotlib import pyplot as plt
import Models2


wr = Models2.MKVOR()

wr2 = Models2.myMod()
frange = np.linspace(0, 1, 100)
#
# plt.plot(frange, [wr.C.eta_s(f) for f in frange])
# plt.plot(frange, [wr.C.eta_o(f) for f in frange])
# plt.plot(frange, [wr.C.eta_r(f) for f in frange])
#
# plt.plot(frange, [wr2.C.eta_s(f) for f in frange])
# plt.plot(frange, [wr2.C.eta_o(f) for f in frange])
# plt.plot(frange, [wr2.C.eta_r(f) for f in frange])
#
# plt.show()
# exit()
wr.dumpAll(hyper=0)

wr.hyper_phi.dumpEos()
wr.hyper_phi.dumpMassesCrust()

wr.hyper_phi_sigma.dumpEos()
wr.hyper_phi_sigma.dumpMassesCrust()

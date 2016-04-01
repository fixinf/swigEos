import eosWrap as eos
import numpy as np
from matplotlib import pyplot as plt
import Models2
from scipy.optimize.minpack import curve_fit


wr = Models2.MKVOR2_fom(0.6)

wr2 = Models2.MKVOR()
frange = np.linspace(0, 1, 1000)
#

wr3 = Models2.myModExp()

# wr.C.fcut_om = 0.6
# wr.C.acut_om = 1.
# wr.C.bcut_om = 100.
#
# wr.C.fcut_rho = 0.67248809
# wr.C.bcut_rho = 23.96047632

def func(x, f, b):
    # print(x, f, b)
    wr.C.fcut_rho = f
    wr.C.bcut_rho = b
    return [wr.C.eta_r(z) for z in x]


# res = curve_fit(func, frange, [wr3.C.eta_r(f) for f in frange])
# print(res)
# exit()

plt.plot(frange, [wr3.C.eta_r(f) for f in frange])

plt.plot(frange, [wr.C.eta_s(f) for f in frange])
plt.plot(frange, [wr.C.eta_o(f) for f in frange])
plt.plot(frange, [wr.C.eta_r(f) for f in frange])

plt.plot(frange, [wr2.C.eta_s(f) for f in frange])
plt.plot(frange, [wr2.C.eta_o(f) for f in frange])
plt.plot(frange, [wr2.C.eta_r(f) for f in frange])

plt.show()
# exit()
wr.dumpAll(hyper=0)

wr.hyper_phi.dumpEos()
wr.hyper_phi.dumpMassesCrust()

wr.hyper_phi_sigma.dumpEos()
wr.hyper_phi_sigma.dumpMassesCrust()

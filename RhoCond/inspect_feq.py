import Models2
from matplotlib import pyplot as plt
import numpy as np

wr = Models2.waleckaMatsui()

m = wr.rcond_nucl

frange = np.linspace(0, 1, 100)
chi = (1 - frange) / np.array([np.sqrt(m.C.eta_r(f)) for f in frange])
# plt.plot(frange, chi)
# plt.show()

# m.loadEos()
# plt.plot(m.nrange/m.n0, m.mu_e)
# plt.plot(m.nrange/m.n0, m.C.m_rho * (1 - m.rho[:, 0]))
# plt.show()
# exit()

m.inspect_f()
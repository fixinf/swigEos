import eosWrap as eos
import numpy as np
from matplotlib import pyplot as plt

# x = np.linspace(0, 6, 100)
# plt.plot(x, np.sin(x))
# plt.show()
# f = 0.65

import Models2

f = 0.65

wr = Models2.KVOR()
# wr.hyper_phi.dumpEos()
# wr.hyper_phi.dumpMassesCrust()
# exit()
m = wr.rcond_nucl
# m.filenames['eos'] += '+'
# m.dumpEos()
m.loadEos()
p_deriv = m.nrange * np.gradient(m._E, m.nrange[1]-m.nrange[0]) - m._E

plt.plot(m.nrange/m.n0, m.rho[:, 0])
plt.plot(m.nrange/m.n0, np.gradient(m.rho[:,0], m.nrange[1]-m.nrange[0]))
plt.show()
# plt.plot(m.nrange/m.n0, m._E)
# plt.plot(m.nrange/m.n0, np.gradient(m._E, m.nrange[1]-m.nrange[0]))
# plt.plot(m.nrange/m.n0, m.nc)
# plt.show()
plt.plot(m.nrange/m.n0, p_deriv*m.mpi4_2_mevfm3)
plt.plot(m.nrange/m.n0, m._P)
plt.show()
exit()
# m.dumpEos()

lc = m.lepton_concentrations().transpose()
nc = m.nc / m.nrange

lines = plt.plot(m.nrange/m.n0, m.concentrations())
plt.legend(lines, m.part_names)
plt.plot(m.nrange/m.n0, nc)
plt.plot(m.nrange/m.n0, lc)

print((m.lepton_concentrations()[0] * m.nrange)[-1])
print(m.nc[-1])

print(m.concentrations()[:, 1] + m.concentrations()[:, 5] - lc[:, 0] - lc[:, 1] - nc)

plt.show()
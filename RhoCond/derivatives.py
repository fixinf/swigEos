import eosWrap as eos
import numpy as np
from matplotlib import pyplot as plt

import Models2

wr = Models2.KVOR()

m = wr.rcond_nucl
m.reset()

mus = m.mu()

plt.plot(m.nrange/m.n0, m.mu_e)
plt.plot(m.nrange/m.n0, np.gradient(m.mu_e, m.nrange[1]-m.nrange[0]))
plt.show()

plt.plot(m.nrange/m.n0, mus)
plt.plot(m.nrange/m.n0, np.gradient(mus[:, 0], m.nrange[1]-m.nrange[0]))
plt.plot(m.nrange/m.n0, np.gradient(mus[:, 1], m.nrange[1]-m.nrange[0]))
plt.show()

plt.plot(m.nrange/m.n0, m._E)
plt.plot(m.nrange/m.n0, np.gradient(m._E, m.nrange[1]-m.nrange[0]))
plt.show()

plt.plot(m.nrange/m.n0, m.Efull())
plt.plot(m.nrange/m.n0, np.gradient(m.Efull(), m.nrange[1]-m.nrange[0]))
plt.show()
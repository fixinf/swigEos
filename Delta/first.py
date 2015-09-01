from scipy.misc.common import derivative
from Models2 import KVOR, _KVOR
import Models2

__author__ = 'const'
import eosWrap as eos
import numpy as np
from matplotlib import pyplot as plt
import Wrapper2

np.set_printoptions(threshold=np.nan)

# C = eos.SCDelta
# C = eos.Walecka
# C = _KVOR
# C = eos.Walecka_d
C = eos.KVOR_d
wr = Wrapper2.Model(C)

# wr = Models2.KVOR()
# wr = Models2.myMod()

m = wr.delta
m.C.SetHyperConstants(2)
n_in = np.array([0.2, 0.25, 0.25, 0, 0, 0, 0, 0, 0, 0.2])
print(eos.mu(n_in, 8, m.C))
print(derivative(lambda z: eos._E(np.array([0.2, 0.25, 0.25, 0, 0, 0, 0, 0, z, .2]), m.C), 0., dx=1e-3))
# exit()

# m.C.set_xs(np.array([1, 1, -28, 30., 30, 30, -18, -18]))
m.verbose = 1
print([i for i in m.C.X_s])
print([i for i in m.C.Q])
print([i for i in m.C.M])
print([i for i in m.C.T])

# plt.plot(wr.sym.nrange, wr.sym.Ebind(wr.sym.nrange))
# plt.show()

m.reset()
print(m.rho)

mu = m.mu()
lines = plt.plot(m.nrange/m.n0, m.concentrations())
line_f = plt.plot(m.nrange/m.n0, m.rho[:,0])

# lines = plt.semilogy(m.nrange*.16/m.n0, m.concentrations())
# plt.ylim([1e-3, 1])
plt.legend(lines + line_f, m.part_names + ['f'], loc=0)
plt.show()

print(mu[:,0] - (mu[:,1] + mu[:,8]))

lines = plt.plot(m.nrange/wr.n0, mu)
plt.legend(lines, m.part_names)
plt.show()
exit()
plt.plot(m.nrange, mu[:,0] + m.mu_e)
Xm, = plt.plot(m.nrange, mu[:,6])
Dm, = plt.plot(m.nrange, mu[:,8])
plt.legend([Xm, Dm], ['Xm', 'Dm'])

plt.show()
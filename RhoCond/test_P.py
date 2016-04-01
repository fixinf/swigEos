import eosWrap as eos
from matplotlib import pyplot as plt
import numpy as np
from scipy.misc.common import derivative
import Models2

# wr = Models2.myModExpOmega(0.66)
wr = Models2.KVOR()
m = wr.rcond_hyper_phi
m2 = wr.nucl
# m2.dumpEos()
m2.loadEos()

m.reset()
plt.plot(wr.nrange/wr.n0, m._P, label='rho')
plt.plot(wr.nrange/wr.n0, m2._P, label='norho')
plt.legend()
plt.show()

pdiff = m.nrange*np.gradient(m._E, m.nrange[1]-m.nrange[0]) - m._E
pdiff *= m.mpi4_2_mevfm3

plt.plot(wr.nrange/wr.n0, m._P, label='rho')
plt.plot(wr.nrange/wr.n0, pdiff, label='rhodiff')
plt.legend()
plt.show()

lines = plt.plot(m.nrange/m.n0, m.concentrations())
plt.legend(lines, m.part_names)
plt.show()

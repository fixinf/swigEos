__author__ = 'const'
import Models2
import eosWrap as eos
import numpy as np
import matplotlib.pyplot as plt

m = Models2.waleckaMatsui()
m2 =  Models2.KVOR()
C = m.sym.C
C2 = m2.sym.C

frange = np.linspace(0, 1, 100)

E = list(map(lambda z: (eos._E(np.array([z, m.n0/2, m.n0/2]), C) / m.n0 - C.M[0]) * 135., frange))
E2 = list(map(lambda z: (eos._E(np.array([z, m.n0/2, m.n0/2]), C2) / m.n0 - C.M[0]) * 135., frange))
plt.plot(1-frange, E)
plt.plot(1-frange, E2)
plt.show()
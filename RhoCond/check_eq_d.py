import eosWrap as eos
import Models2
import numpy as np
from matplotlib import pyplot as plt

wr = Models2.MKVOR_d()

m = wr.delta_sym
m.loadEos()

i = 150
r = m.rho[i, :]
n =m.nrange[i]
print(r)
print(m.nrange[i], 4*r[10] + 2*r[1])

print(m.eq([4*r[10]], m.nrange[i], r[0]))
plt.plot(m.nrange/m.n0, [m.eq([z], m.nrange[i], r[0]) for z in m.nrange])
plt.ylim([-5., 5.])
plt.show()
m.inspect_eq()
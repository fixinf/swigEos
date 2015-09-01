import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join

C = Models.Cubero()
wr = Wrapper(C)

C2 = Models.Cubero_cut03()
wr2 = Wrapper(C2)

print(wr.m_pi * C.M[0])
n = np.linspace(0, 3, 100)
print(eos.f_eq(np.array([wr.n0/2, wr.n0/2])/1, np.array([0.]), 1, C))
print(wr.K())
E, f = wr.Esymm(n, ret_f=1)
E2, f2 = wr2.Esymm(n, ret_f=1)
frange = np.linspace(0, 1, 100)
plt.plot(frange, [[C.U(z), C2.U(z)] for z in frange])
plt.show()
plt.plot(n/wr.n0, np.array([f, f2]).transpose())
plt.show()
plt.plot(n/wr.n0, (E/n - C.M[0])*wr.m_pi, n/wr.n0, wr2.ESbind(n))
plt.show()
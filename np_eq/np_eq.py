import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join
from pylab import pi

def np_eq(n, J, L, n0):
    eps = (n - n0)/n0
    return (n/6.) * (1 - (1 + 128 * (J + L*eps/3.)**(3) / (pi**2 * n) )**(-1))

C = Models.KVOR()
wr = Wrapper(C)



J = wr.J()/wr.m_pi
L = wr.L()/wr.m_pi

print J*wr.m_pi, wr.Jdiff(wr.n0)
print wr.L(), wr.Ldiff()
print wr.Ksymm(), wr.Ksymm_diff()
exit()

wr.reset()

rho = wr.concentrations()

np_eq = map(lambda z: np_eq(z, J, L, wr.n0), wr.n)

plt.plot(wr.n/wr.n0, rho[:, 1], wr.n/wr.n0, np_eq)
plt.show()
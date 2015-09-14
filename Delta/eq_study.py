__author__ = 'const'
from scipy.misc.common import derivative
from Models2 import KVOR, _KVOR
import Models2

__author__ = 'const'
import eosWrap as eos
import numpy as np
from matplotlib import pyplot as plt
import Wrapper2
from numpy import array as arr

np.set_printoptions(threshold=np.nan)

wr = Models2.MKVOR_d()

m = wr.delta_phi
m.inspect_f()
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
exit()
#void func_f_eq(double * p, double * hx, int m, int _n, void * adata)
out = eos.dArray(1)
n = eos.dArray(8)
_n = [1.13689708e+00,   3.34609432e-01,   2.77758946e-01,
    5.30624945e-16,   1.97327183e-16,   0.00000000e+00,   1.53208617e-01,
    0.00000000e+00,   5.13019741e-02]
21
for i, __n in enumerate(_n):
    n[i] = __n

params = eos.func_f_eq_params()
params.n = n
params.C = m.C
params.dimN = 8
params.df = 1e-3

def set_arg(x):
    a = eos.dArray(1)
    a[0] = x
    return a

def func(f):
    eos.func_f_eq(set_arg(f), out, 1, 1, params)
    return out[0]

frange = np.linspace(0, 1., 100)
plt.plot(frange, list(map(func, frange)))
plt.show()
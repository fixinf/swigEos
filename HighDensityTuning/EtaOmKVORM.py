import matplotlib
from scipy.misc.common import derivative
from time import sleep
from math import exp
matplotlib.use('QT4Agg')
from matplotlib.widgets import Slider
import numpy as np
import matplotlib.pyplot as plt
import Models
from Wrapper import Wrapper
from scipy import optimize

C = Models.myMod()
C1 = Models.myMod()

wr = Wrapper(C)
wr1 = Wrapper(C1)

C.omega_c = 0

C1.omega_c = 0

C1.omega_kind = 2
C1.omega_a = 0.15669423
C1.omega_b = 7.09984178
C1.omega_f = 0.95

npoints = 10
f_solve = np.linspace(0.5, 1., npoints)

def func(x):
    C1.omega_a = x[0]
    C1.omega_b = x[1]
    res = map(lambda z: C.eta_o(z) - C1.eta_o(z), f_solve)
    res = np.array(res)
    return np.sum(res**2)

# res = optimize.minimize(func, [C1.omega_a, C1.omega_b])
# print res.x

frange = np.linspace(0, 1, 100)

plt.plot(frange, map(C.eta_o, frange), frange, map(C1.eta_o, frange))
plt.ylim([0.5, 1.5])
plt.show()

n = np.linspace(0, 8*wr.n0, 100)

Es, fs = wr.Esymm(n, ret_f=1)
plt.plot(n/wr.n0, fs)
plt.show()

wr1.testDanielewicz()
# wr1.testPodsiedlowski()
wr1.dumpVs()

wr.reset()
wr1.reset()

wr.setDriver()
wr1.setDriver()


n,m,r,mg= wr.stars()
n1,m1,r1,mg1= wr1.stars()

print max(m), max(m1)



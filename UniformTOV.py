import matplotlib
from scipy.misc.common import derivative
from time import sleep
from math import exp, pi
from cmath import sqrt
from scipy.interpolate.interpolate import interp1d
from numpy import arcsin
matplotlib.use('QT4Agg')
from matplotlib.widgets import Slider
import numpy as np
import matplotlib.pyplot as plt
import Models
from Wrapper import Wrapper
from scipy import optimize


C = Models.KVOR()
wr = Wrapper(C)
wr.reset()
wr.setDriver()
n, m, r, mg1 = wr.stars(nmax=3.9)
iE = interp1d(wr.n, wr.E)
iP = interp1d(wr.n, wr.P)
Ec = iE(n)
Pc = iP(n)
res = []
res2 = []

#Test M/R ratio 
for i, _n in enumerate(n):
    Emean = m[i] * 1.4766/ (4 * pi * r[i]**3 / 3) # km^{-2}
    Emean /= 2.6115e-4 # fm^{-4}www.win2.n/
    print Ec[i], Emean 
    res.append(2 * m[i] * 1.4766/r[i])
    res2.append(1 - ((Pc[i] + Emean)/(3 * Pc[i] + Emean))**2)

res = np.array(res)
res2 = np.array(res2)
plt.plot(n/wr.n0, res)
plt.plot(n/wr.n0, res2, ls='--')
plt.plot(n/wr.n0, (res2-res)/res)
plt.show()

#Test baryon number
C.Hyper = 0
n, m, r, mg1, mg2 = wr.stars_crust()
plt.plot(n/wr.n0, mg2)

A = 2 * m * 1.4766 / r
N_B = n * r**3 / (2 * A**(1.5)) * (arcsin(np.sqrt(A)) - np.sqrt(A*(1-A)))
plt.plot(n/wr.n0, (0.0004898007281478712)*N_B)

plt.show()

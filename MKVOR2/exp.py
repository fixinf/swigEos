from cmath import cosh

from scipy.optimize.minpack import leastsq

import Models2
import Models
from matplotlib import pyplot as plt
import numpy as np

wr1 = Models2.myMod()
wr2 = Models2.myMod()

m1 = wr1.nucl
m2 = wr2.nucl

C1 = m1.C
C = Models.myModR()

C.SetHyperConstants(2)

C.Csp = 380.

C.phi_gamma = 3.0
C.phi_z=3.5

frange = np.linspace(0, 1, 100)
# plt.plot(frange, [C1.eta_r(f) for f in frange])

def func(f, x):
    beta, gamma, rho_f, rho_a = x
    res = ((1 + beta*f**2) / (1 + beta*C.f0**2))**gamma
    if f > rho_f:
        res /= cosh(rho_a * (f - rho_f)**2)
    return res.real

frange_to_fit = frange[:60]

def f_to_fit(x):
    return np.array([func(f, x) - C1.eta_r(f) for f in frange_to_fit])
x0 = [C.beta, C.gamma, C.rho_f, C.rho_a]
print(f_to_fit(x0))
out = leastsq(f_to_fit, [C.beta, C.gamma, C.rho_f, C.rho_a])
print(out)
plt.plot(frange, [C1.eta_r(f) for f in frange])
# plt.plot(frange, [func(f, out[0]) for f in frange])
C3 = Models2._myModExp()
plt.plot(frange, [C3.eta_r(f) for f in frange])
plt.show()

C3 = Models2._myModExp()
wr3 = Models2.myModExp()
wr3.hyper.reset(timeout=6)
plt.plot(wr3.nrange[:wr3.hyper.rho.shape[0]]/wr3.hyper.n0, wr3.hyper.concentrations())
plt.show()



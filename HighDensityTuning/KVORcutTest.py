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

model = Models.KVOR_cut_02
C1 = Models.KVOR_tan_02()
C = model()
wr = Wrapper(C)
wr1 = Wrapper(C1)

# C.omega_kind = 2
# C.omega_b = 100
# C.omega_a = -0.5
# C.omega_f += 0.05

# def f_cut(x):
#     C.omega_b = x[0]
#     C.omega_f = x[1]
#     res = [C.eta_o(C.omega_f) - C1.eta_o(C.omega_f),
#             derivative(C.eta_o, C.omega_f, dx=1e-3) - derivative(C1.eta_o, C.omega_f, dx=1e-3)]
#     res = np.array(res)
#     return np.sum(res*res)
# 
# print f_cut([C.omega_b, C.omega_f])
# 
# res = optimize.minimize(f_cut, [C.omega_b, C.omega_f])
# 
# print res.x

frange = np.linspace(0, 1, 100)
plt.plot(frange, map(C.eta_o, frange), frange, map(C1.eta_o, frange))
plt.show()

C.phi_meson = 0
wr.reset(hyper=1, npoints=100)
wr1.reset(hyper=1, npoints=100)

plt.plot(wr.n/wr.n0, wr.concentrations(), wr.n/wr.n0, wr1.concentrations())
plt.show()
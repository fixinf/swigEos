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

C1.omega_c = 0

C1.rho_tan_b = 20
C1.rho_tan_a = 3
C1.rho_sat_f1 = 0.3

frange = np.linspace(0, 1, 100)

plt.plot(frange, list(map(C.eta_r, frange)), frange, list(map(C1.eta_r, frange)))
# plt.ylim([0, 50])
plt.show()

wr.reset()
wr1.reset()

wr.setDriver()
wr1.setDriver()

n,m,r,mg= wr.stars()
n1,m1,r1,mg1= wr1.stars()

print(max(m), max(m1))



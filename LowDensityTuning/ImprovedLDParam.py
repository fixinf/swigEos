import matplotlib
from scipy.misc.common import derivative
from time import sleep
from math import exp
from scipy.interpolate.interpolate import interp1d
from matplotlib.widgets import Slider
matplotlib.use('QT4Agg')
import numpy as np
import matplotlib.pyplot as plt
import Models
from Wrapper import Wrapper
from scipy import optimize
import eosWrap as eos

C1 = Models.myMod()
wr1 = Wrapper(C1);
n = np.linspace(0, wr1.n0, 100)

data = np.loadtxt('/home/const/GrabbedFigures/HebelerSchwenk2014/HebelerSchwenk_Out.dat',
                   skiprows=1)

C = eos.ImprovedLDParam();
wrLD= Wrapper(C)

C.Cs = 234.1580555799
C.Co = 134.8826104616
C.b = 0.0046776700
C.c = -0.0029781609
C.Cr = 81.7485399416

C.alpha = 0.4
C.z = 0.65

C.omega_a = 0.8
C.omega_f = 0.55

C.beta = 1.2
C.gamma = 7.5
          
C.phi_a = -0.0
C.phi_f = 0.28

C.d = -0.5

C.rho_f = 0.75
C.rho_a = 1000

C.f0 = 0.27
C.rho_kind = 1
C.rho_power = 2.0
C.omega_c = -20000

C.SetHyperConstants(2)

C.Csp = 380.

C.phi_gamma = 3.0
C.phi_z=3.5

C.rho_f = 0.45
C.rho_a = 100.

frange = np.linspace(0, 1, 100)


# wrLD.solve(f0=0.27, K0=240., J0=30.)

C.ar = 1.94841242712
C.br= 1.53538400433
C.cr = 1.91620264189

# C.ar, C.br, C.cr = [1.9623089869889265, 1.4986992352213222, 1.9389911745612569]
C.ar, C.br, C.cr = [2.0685962608175834, 1.1117463744674665, 2.2196563820735964]

C.f_stop = 0.237591836735
C.f_stop_denom = C.f_stop
C.a_sigma = -0.0
C.power_0 = 0
C.power_1 = 3

plt.plot(frange, map(C.eta_r, frange), frange, map(C1.eta_r, frange))
plt.plot(frange, map(C.eta_o, frange), frange, map(C1.eta_o, frange))
plt.plot(frange, map(C.eta_s, frange), frange, map(C1.eta_s, frange))
plt.ylim([0, 5])
print C.eta_r(0.23)
plt.show()

plt.plot(data[:,0], data[:, 1:],c='black',lw='1')
plt.plot(n/wr1.n0, 135*(wrLD.Eneutr(n)/n - C1.M[0]),c='blue')
plt.plot(n/wr1.n0, 135*(wr1.Eneutr(n)/n - C1.M[0]),c='blue')
plt.show()

wrLD.dumpVs()
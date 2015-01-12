import matplotlib
from scipy.misc.common import derivative
from time import sleep
from math import exp
from scipy.interpolate.interpolate import interp1d
from scipy import interpolate
from matplotlib.widgets import Slider
from numpy import sqrt
matplotlib.use('QT4Agg')
import numpy as np
import matplotlib.pyplot as plt
import Models
from Wrapper import Wrapper
from scipy import optimize
import eosWrap as eos

C = eos.InterpolatedScalings()
C.debug = 0
C1 = Models.myMod()
wr = Wrapper(C)
wr1 = Wrapper(C1)

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

npoints = 100
mix = 0.9



a_s = -0.0

frange = np.linspace(-1e-2, 1+1e-3, npoints)
etar = np.array(map(C1.eta_r, frange))
etao = np.array(map(C1.eta_o, frange))
u = np.array(map(C1.U, frange))
etas = np.array(map(C1.eta_s, frange))


n_sparse = 30
f_sparse = np.linspace(-1e-2, 1+1e-3, n_sparse)

etar_sparse = np.array(map(C1.eta_r, f_sparse)) 

start = n_sparse*0.56
stop = n_sparse

print  [C1.eta_r(f_sparse[start]) for i in f_sparse[start:stop]]
etar_sparse[start:stop] = [C1.eta_r(f_sparse[start]) - sqrt(i) for i in f_sparse[start:stop]]

plt.plot(f_sparse, etar_sparse)
plt.show()

print f_sparse.shape, etar_sparse.shape

# i_etar = interp1d(f_sparse, etar_sparse, kind='quadratic')
i_etar = interpolate.UnivariateSpline(f_sparse, etar_sparse, s=1)
etar_dense = i_etar(frange)

plt.plot(frange, etar_dense)
plt.show()

C.rho_akima = 0

C.set_eta_o(frange, etao)
C.set_eta_r(frange, etar_dense)
C.set_eta_s(frange, etas)
C.set_U(frange, u)

plt.plot(frange, map(C.eta_r, frange), frange, map(C1.eta_r, frange))
plt.plot(frange, map(C.eta_o, frange), frange, map(C1.eta_o, frange))
plt.plot(frange, map(C.eta_s, frange), frange, map(C1.eta_s, frange))
plt.ylim([0, 5])
print C.eta_r(0.23)
plt.show()

wr.dumpVs()
# wr.set=1
wr.setDriver()
n, m, r, mg = wr.stars()
print 'M_NS_MAX =', max(m)
sleep(3)

wr1.C.phi_meson=0
wr1.reset(npoints=npoints, timeout=2, hyper=1)



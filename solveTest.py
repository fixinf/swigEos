import eosWrap as eos
from Wrapper import Wrapper
import matplotlib
matplotlib.use("QT4Agg")
from pylab import *
from numpy import array, linspace
from matplotlib.widgets import Slider, CheckButtons, RadioButtons

# x = linspace(0.0, 3.0, 100)
# l = plot(x, sin(x),x, sin(2*x),x, sin(3*x))
# colors = [c.properties()['color'] for c in l]
# print colors
# show()

C = eos.KVOR_mod()
C.SetHyperConstants(2)
# 
C.alpha = 0.85
C.z = 0.65
   
C.omega_a = 1.67
C.omega_f = 0.44
   
C.rho_a = 0.0
   
C.phi_a = -1.01
C.phi_f = 0.33
   
C.d = -5.0
   
C.phi_gamma = 3.0
C.phi_z = 3.5
   
C.sprime = 0
C.Csp = 380.0
   
# C.beta = 4.54
# C.gamma = 3.78

f0 = 0.26


print C.X_s[7]

print C.X_s[7]


wr = Wrapper(C)
C.Csp = 1.0
n0 = wr.n0
K0 = 275.0
print eos.EBind(np.array([0.195, wr.n0/2, wr.n0/2]), C)
print eos.f_eq(np.array([n0/2, n0/2]), np.array([0.0]), 1, C)
wr.solve()
wr.solve(f0=0.26)
pause(10)
for f in linspace(C.f0, f0, 20):
    C.f0 = f
    wr.solve(f0=C.f0, K0=K0) 
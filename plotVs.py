import eosWrap as eos
from Wrapper import Wrapper
from pylab import *
from numpy import array, linspace
from matplotlib.widgets import Slider, CheckButtons, RadioButtons
from scipy.misc.common import derivative
# x = linspace(0.0, 3.0, 100)
# l = plot(x, sin(x),x, sin(2*x),x, sin(3*x))
# colors = [c.properties()['color'] for c in l]
# print colors
# show()

C = eos.KVOR_mod()
C.SetHyperConstants(2)
C.Csp = 1.0
C.alpha = 0.85
C.omega_a = 12.0
C.rho_a = 0.0
C.phi_f = 0.35
C.phi_a = -3.5
C.omega_f = 0.55
C.d = -9.8
C.phi_gamma = -0.8
C.z = 0.65
C.phi_z = 3.1

C.beta = 4.54
C.gamma = 3.78

C2 = eos.KVOR_mod()
C2.Csp = 1.0
C2.SetHyperConstants(2)
C2.alpha = 0.85
C2.omega_a = 1.67
C2.rho_a = 0.0
C2.phi_f = 0.33
C2.phi_a = -1.01
C2.omega_f = 0.44
C2.d = -4.21
C2.phi_gamma = -0.8
C2.z = 0.65
C2.phi_z = 3.1

C2.beta = 4.54
C2.gamma = 3.78

f0 = 0.26

wr = Wrapper(C)
wr2 = Wrapper(C2)

n_p = linspace(0.0, 4.0, 1000)

wrlist = [wr, wr2]

for wr in wrlist:
    for f in linspace(0.195, f0, 200):
        wr.C.f0 = f
        wr.solve(f0=f) 
        
    print 'Done!'
    pause(1)

print eos.J(wr.n0, C2)
print 3*derivative(lambda z: eos.J(z, C2), wr.n0, dx=1e-3)
pause(5)

fig, ax = plt.subplots(1,2)

UpperY = [15.639, 60, 140, 210] 
UpperX = [2, 2.75, 3.5, 4.5]
LowerX= [2, 2.5, 3.5, 4, 4.5]
LowerY = [7.5, 17.49, 40.0, 60.0, 60.0]

ax[0].plot(UpperX, UpperY, c='red')
ax[0].plot(LowerX, LowerY, c='red')
clist = ['blue', 'green']
for i,wr in enumerate(wrlist):
    _P = wr.Psymm(n_p)
    ax[0].semilogy(n_p[1:]/wr.n0, _P, c=clist[i])
    _E = wr.Esymm(n_p)[1:]
    Vs = np.diff(_P)/(np.diff(_E)*wr.m_pi**4*wr.const)
    ax[1].plot(n_p[2:]/wr.n0, Vs)


ax[0].set_xlabel(r'$n/n_0$')
ax[0].set_ylabel(r'$P_{symm}, \quad MeV/fm^{3}$')

ax[1].set_xlabel(r'$n/n_0$')
ax[1].set_ylabel(r'$v_{s}^2, \, SNM$')

    
show()
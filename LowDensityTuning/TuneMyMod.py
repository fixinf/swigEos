import matplotlib
matplotlib.use('Qt4Agg')
from Wrapper import Wrapper
from matplotlib.widgets import Slider
import numpy as np
import Models
import eosWrap as eos
import matplotlib.pyplot as plt

C = Models.myModLow()
wr = Wrapper(C)

frange = np.linspace(0., 1., 500)

C.omega_f_low = 0.2
C.omega_a_low = 9.44

C.Delta = 0*25.23/135./wr.n0

C.rho_f_low = 0.2
C.rho_a_low = 3.

n = np.linspace(0, 1., 100)
E = wr.Eneutr(n)
# wr.solve(f0=C.f0, J0=32., K0=240.)
x = []
low = []
up = []

with open('HebelerSchwenk_Out.dat', 'r') as f:
    for line in f:
#         print line.strip(' ').split()
        _x, _l, _u = line.split()
        x.append(float(_x))
        low.append(float(_l))
        up.append(float(_u))
fig, ax = plt.subplots(1, 4)
plt.subplots_adjust(left=0.25, bottom = 0.25)

ebindN = (E/n - C.M[0])*135.
l, = ax[0].plot(n/wr.n0, ebindN)

ax[0].plot(x, low, x, up)

ax[0].set_xlim([0., 1.])
ax[0].set_ylim([0., 20])

ax[1].set_ylim([-2, 5])

lom, = ax[1].plot(frange, map(C.eta_o, frange))
lr, = ax[1].plot(frange, map(C.eta_r, frange))

lsymm, = ax[2].plot(n/wr.n0, (wr.Esymm(n)/n - C.M[0])*wr.m_pi)
ax[2].set_xlim([0., 1.])

EN = wr.Eneutr(n)
PN = wr.P_N(n)
lvs,=ax[3].plot(n[1:]/wr.n0, np.diff(PN)/np.diff(EN)/wr.const/wr.m_pi**4)

axRhoF = plt.axes([0.25, 0.1, 0.3, 0.03])
axOmegaF = plt.axes([0.25, 0.15, 0.3, 0.03])
axRhoA = plt.axes([0.6, 0.1, 0.3, 0.03])
axOmegaA = plt.axes([0.6, 0.15, 0.3, 0.03])

slRhoF = Slider(axRhoF, 'Rho F', 0., C.f0, valinit = C.rho_f_low)
slOmegaF = Slider(axOmegaF, 'Omega F', 0., C.f0, valinit = C.omega_f_low)
# slRhoA = Slider(axRhoA, 'Rho A', -1000., 30., valinit = C.rho_a_low)
slRhoA = Slider(axRhoA, 'Rho A', 1.5, 5., valinit = C.rho_a_low)
slOmegaA = Slider(axOmegaA, 'Omega A', -30., 30., valinit = C.omega_a_low)

sliders = [slRhoA, slRhoF, slOmegaA, slOmegaF]

def update(val):
    C.rho_f_low = slRhoF.val
    C.rho_a_low = slRhoA.val
    C.omega_a_low = slOmegaA.val
    C.omega_f_low = slOmegaF.val
    E = (wr.Eneutr(n)/n - C.M[0])*135.
    l.set_ydata(E)
    lr.set_ydata(map(C.eta_r, frange))
    lom.set_ydata(map(C.eta_o, frange))
    lsymm.set_ydata((wr.Esymm(n)/n - C.M[0])*wr.m_pi)
    EN = wr.Eneutr(n)
    PN = wr.P_N(n)
    lvs.set_ydata(np.diff(PN)/np.diff(EN)/wr.const/wr.m_pi**4)

    
for slider in sliders:
    slider.on_changed(update)

plt.show()
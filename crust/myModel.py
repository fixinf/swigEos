import matplotlib
matplotlib.use('QT4Agg')
import Models
import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt

C = Models.KVOR()
wr = Wrapper(C)
N, M, R = wr.stars_crust(ncut_crust=0.5, inter='linear', ncut_eos=0.6, nmin=0.4, nmax=4., npoints=50, show=1)

Mk = []
Rk = []

with open('KVOR_BPS_Klahn', 'r') as f:
    for line in f:
        r, m = line.split()
        Rk.append(float(r))
        Mk.append(float(m))
      
plt.plot(Rk, Mk)
plt.plot(R, M)
plt.xlabel(r'$R \, [km]$', fontsize=18)
plt.ylabel(r'$M/M_\odot$', fontsize=18)
plt.show()

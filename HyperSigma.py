import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join
from os import path
import os

model = Models.KVOR_cut_03
C = model()
C.set_xs(np.array([0., 0., -28., 30., 30., 30., -15., -15.]))
wr = Wrapper(C)

C.sigma_kind = 1

zl = 3.5
zx = 3.5

alpha = 5.5

# wr.solve(f0=C.f0, J0=28., K0=240.)

folder = '/home/const/Dropbox/Documents/For DN/XLambda/data/hyper_alpha=%i_%s_%.2i_%.2i'%(alpha,
                                                                          model.__name__,
                                                                  zl,
                                                                  zx)
if not path.exists(folder):
    os.makedirs(folder)

C.set_hs_alpha(np.array([0., 0.] + [alpha for i in xrange(6)]))
C.set_hs_z(np.array([0., 0., zl, zl, zl, zl, zl, zl]))

print 'hey'

frange = np.linspace(0., 1., 100)
scF = map(lambda z: ((1 + zl * C.f0) / (1 + zl * z))**alpha, frange)
plt.plot(frange, scF)
plt.show()

wr.reset(hyper=1, nmin=0., nmax=8*wr.n0, npoints=200, timeout=3)
wr.setDriver()

plt.plot(wr.n/wr.n0, wr.concentrations(), wr.n/wr.n0, wr.rho[:,0])
plt.show()

scN = map(lambda z: ((1 + zl * C.f0) / (1 + zl * z))**alpha, wr.rho[:, 0])

tabscF = np.array([frange, scF]).transpose()
tabscN = np.array([wr.n/wr.n0, scN]).transpose()

tabscF = tabulate(tabscF, ['f', 'x_H'], tablefmt='plain')
tabscN = tabulate(tabscN, ['f', 'x_H'], tablefmt='plain')

with open(join(folder, 'scF.dat'), 'w') as f:
    f.write(tabscF)

with open(join(folder, 'scN.dat'), 'w') as f:
    f.write(tabscN)

n, m, r, mg1 = wr.stars(npoints=30)

plt.plot(n/wr.n0, m)
plt.show()

mtable = tabulate(np.array([n/wr.n0, m, r]).transpose(), ['n/n_0', 'M/M_sun', 'R [km]'], tablefmt='plain')
with open(join(folder, 'masses.dat'), 'w') as f:
    f.write(mtable)

rho = wr.concentrations()
print wr.n.shape
print rho.shape
rtab = []
for i, _n in enumerate(wr.n/wr.n0):
#     print _n
#     print rho[i]
    rtab.append(np.concatenate(([_n], rho[i])))

rtable = tabulate(rtab, ['n/n_0', 'n', 'p', 'L', 'S-', 'S0', 'S+', 'X-', 'X+'],
                  tablefmt='plain')

with open(join(folder, 'rho.dat'), 'w') as f:
    f.write(rtable)

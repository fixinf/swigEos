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

zl = 1.
zx = 1.

alpha = 0.

# wr.solve(f0=C.f0, J0=28., K0=240.)

C.set_hs_alpha(np.array([0., 0.] + [alpha for i in xrange(6)]))
C.set_hs_z(np.array([0., 0., zl, 0., 0., 0., zx, zx]))

print 'hey'
C.phi_meson = 1
# wr.reset(hyper=1, nmin=0., nmax=8.5*wr.n0, npoints=200)

C.phi_gamma = 4.5
C.phi_z = 5.5

wr.setDriver()
n, m, r, mg1, mg2 = wr.stars_crust_hyper(npoints=30)
print max(m)
# exit()
rhos = wr.concentrations()
# print rhos
plt.plot(wr.n/wr.n0, rhos)
plt.show()

exit()
frange = np.linspace(0., 1., 100)

etapNS = map(C.eta_p, wr.rho[:,0])

tabetap = tabulate(np.array([wr.n/wr.n0, etapNS]).transpose(), 
                   ['n', 'eta_phi'], tablefmt='plain')


scF = map(lambda z: ((1 + zl * C.f0) / (1 + zl * z))**alpha, frange)
scN = map(lambda z: ((1 + zl * C.f0) / (1 + zl * z))**alpha, wr.rho[:, 0])

fig, ax = plt.subplots(2, 1)
ax[0].plot(frange, scF, frange, map(C.eta_p, frange))
ax[1].plot(wr.n/wr.n0, etapNS, wr.n/wr.n0, map(C.eta_p, wr.rho[:,0]))
plt.show()
tabscF = np.array([frange, scF]).transpose()
tabscN = np.array([wr.n/wr.n0, scN]).transpose()

tabscF = tabulate(tabscF, ['f', 'x_H'], tablefmt='plain')
tabscN = tabulate(tabscN, ['f', 'x_H'], tablefmt='plain')

folder = '/home/const/Dropbox/Documents/For DN/XLambdaPhi/data/hyper_alpha=%.1f_%s_%.1f_%.1f_%.1f_%.1f'%(alpha,
                                                                          model.__name__,
                                                                  zl,
                                                                  zx,
                                                                  C.phi_z,
                                                                  C.phi_gamma)
# with open(join(folder, 'eta_p_n.dat'), 'w') as f:
#     f.write(tabetap)

if not path.exists(folder):
    os.makedirs(folder)



# with open(join(folder, 'scF.dat'), 'w') as f:
#     f.write(tabscF)
# 
# with open(join(folder, 'scN.dat'), 'w') as f:
#     f.write(tabscN)

# wr.dumpScalings(folderName=folder)

# wr.dumpHyper(folderName=folder)

# mtable = tabulate(np.array([n/wr.n0, m, r]).transpose(), ['n/n_0', 'M/M_sun', 'R [km]'], tablefmt='plain')
# with open(join(folder, 'masses.dat'), 'w') as f:
#     f.write(mtable)

# rho = wr.concentrations()
# print wr.n.shape
# print rho.shape
# rtab = []
# for i, _n in enumerate(wr.n/wr.n0):
# #     print _n
# #     print rho[i]
#     rtab.append(np.concatenate(([_n], rho[i])))
# 
# rtable = tabulate(rtab, ['n/n_0', 'n', 'p', 'L', 'S-', 'S0', 'S+', 'X-', 'X+'],
#                   tablefmt='plain')
# 
# with open(join(folder, 'rho.dat'), 'w') as f:
#     f.write(rtable)

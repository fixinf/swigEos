import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join
from os import path
import os

model = Models.KVOR
C = model()
C.set_xs(np.array([0., 0., -28., 30., 30., 30., -15., -15.]))
wr = Wrapper(C)

C.sigma_kind = 1

zl = 2.
zx = 2.

alpha = 2.

folder = '/home/const/Dropbox/Documents/For DN/XLambda/data/alpha=%i_%s_%.2i_%.2i'%(alpha,
                                                                          model.__name__,
                                                                  zl,
                                                                  zx)
if not path.exists(folder):
    os.makedirs(folder)

C.set_hs_alpha(np.array([0., 0.] + [alpha for i in xrange(6)]))
C.set_hs_z(np.array([0., 0., zl, 0., 0., 0., zx, zx]))

print 'hey'

wr.reset(hyper=1, nmin=0., nmax=8*wr.n0, npoints=200)
wr.setDriver()
n, m, r, mg = wr.stars(npoints=50)

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


def func(zl, zx):
    C.set_hs_z(np.array([0., 0., zl, 0., 0., 0., zx, zx]))
    wr.reset(hyper=1, nmin=0., nmax=8*wr.n0, npoints=800)
    wr.setDriver()
    n, m, r, mg = wr.stars(npoints=50)
    fig, ax = plt.subplots(3)
    ax[0].plot(wr.n/wr.n0, wr.concentrations())
    ax[1].plot(wr.n/wr.n0, wr.rho[:, 0])
    lines = ax[1].plot(wr.n/wr.n0, map(lambda z: [C.Xs(2, z), C.Xs(6, z)],
                                        wr.rho[:,0]))
    ax[1].legend(lines, ['L', 'S'], loc=0)
    ax[2].plot(n/wr.n0, m)
    plt.show()
    return np.max(m)
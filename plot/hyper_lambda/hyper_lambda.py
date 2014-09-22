#!/usr/bin/python
import eosWrap as eos
import matplotlib
from scipy.misc.common import derivative
matplotlib.use('QT4Agg')
import matplotlib.pyplot as plt
from Wrapper import Wrapper
import numpy as np
import os
from pylab import pause
from scipy import interpolate
imgfolder = os.path.join('..','img')

C = eos.KVOR_mod()
C.SetHyperConstants(2)
wr = Wrapper(C)
n = np.linspace(0.0, 4.0, 2000)

uX = []
uY = []

def set(params):
    global C
    for par, val in params.iteritems():
        execstring = 'C.'+par+'='+str(val)
        print execstring
        exec(execstring)

# params = dict(
#     alpha=0.85,
#     d = -4.96,
#     f0 = 0.27,
#     a_p = -0.8,
#     f_p = 0.27,
#     a_o = 6.37,
#     f_o = 0.53
# )



C.Csp = 1

# C.Cs = 179.56233875157545
# C.Co =  87.59963973682763
# C.Cr = 100.63642242792484
# C.b = 0.007734608051455927
# C.c = 0.0003446178665624873

C.Cs = 179.56233875157545
C.Co =  87.599595
C.Cr = 100.266691
C.b = 0.008860 
C.c = -0.003312 

C.alpha = 0.85
C.z = 0.65
    
C.omega_a = 6.45
C.omega_f = 0.53
    
C.rho_a = 0.0
    
C.phi_a = -0.85
C.phi_f = 0.28
    
C.d = -5.5
    
C.phi_gamma = 3.0
C.phi_z = 3.5
 
C.sprime = 0
C.Csp = 380.0
    
C.beta = 0.5
C.gamma = 8.0
 

C.f0 = 0.27

C.SetHyperConstants(2)
print [i for i in C.X_s]
C.set_xs(np.array([0.0, 0.0, -30.0, 30.0, 30., 30., -18., -18.]))
print [i for i in C.X_s]

wr.solve(f0=C.f0, iter=3000)

npoints = 500

rhos_lambda = []
masses_lambda = []
nmax = 5.0
C.phi_meson = 1
# ULlist = [-30, -28, -25, -20, -15, -10, -5, 0.0]
ULlist = np.linspace(-30.0, -10.0, 20)
for UL in ULlist:
    rho_crit = [0 for i in xrange(8)]
    change = [1 for i in xrange(8)]
    C.set_xs(np.array([1., 1., UL, 30., 30., 30., -18., -18.]))
    wr.reset(hyper=1, nmax=nmax, npoints=npoints, iter=100)
#     plt.plot(wr.n/wr.n0, wr.rho[:,1:])
#     plt.show()
    wr.setDriver()
    N, M, R = wr.stars(nmax=nmax, npoints=npoints)
    masses_lambda.append(np.max(M))
    for i, rho in enumerate(wr.rho):
#         print rho
#         print rho[1:]
        for j, n in enumerate(rho[1:]):
            if n > 1e-6:
                if change[j]:
                    rho_crit[j] = wr.n[i]
                    change[j] = 0
    rhos_lambda.append(rho_crit)


rhos_xi = []
masses_xi=[]
# UXlist = [-18, -15, -10, -5, 0.0]
UXlist = np.linspace(-18.0, -5.0, 20)
for UX in UXlist:
    rho_crit = [0 for i in xrange(8)]
    change = [1 for i in xrange(8)]
    C.set_xs(np.array([1., 1., -30, 30., 30., 30., UX, UX]))
    wr.reset(hyper=1, nmax=nmax, npoints=npoints, iter=100)
    wr.setDriver()
    N, M, R = wr.stars(nmax=nmax, npoints=npoints)
    masses_xi.append(np.max(M))
    for i, rho in enumerate(wr.rho):
#         print rho
#         print rho[1:]
        for j, n in enumerate(rho[1:]):
            if n > 1e-6:
                if change[j]:
                    rho_crit[j] = wr.n[i]
                    change[j] = 0
    rhos_xi.append(rho_crit)


rhos_lambda = np.array(rhos_lambda)/wr.n0
rhos_xi = np.array(rhos_xi)/wr.n0
print rhos_lambda
print rhos_xi
with open('n_c_lambda_%.2f_%i_%i.dat'%(C.f0, C.phi_meson, C.sprime), 'w') as f:
    for i, U in enumerate(ULlist):
        f.write('%f   '%U)
        for n_c in rhos_lambda[i, 2:]:
            f.write('%f  '%n_c)
        f.write('\n')

with open('n_c_xi_%.2f_%i_%i.dat'%(C.f0, C.phi_meson, C.sprime), 'w') as f:
    for i, U in enumerate(UXlist):
        f.write('%f   '%U)
        for n_c in rhos_xi[i, 2:]:
            f.write('%f  '%n_c)
        f.write('\n')

with open('masses_lambda_%.2f_%i_%i.dat'%(C.f0, C.phi_meson, C.sprime), 'w') as f:
    for i, U in enumerate(ULlist):
        f.write('%f   '%U)
        f.write('%f  '% masses_lambda[i])
        f.write('\n')

with open('masses_xi_%.2f_%i_%i.dat'%(C.f0, C.phi_meson, C.sprime), 'w') as f:
    for i, U in enumerate(UXlist):
        f.write('%f   '%U)
        f.write('%f  '%masses_xi[i])
        f.write('\n')

fig, ax = plt.subplots(1, 2)
linesL = ax[0].plot(ULlist, rhos_lambda[:, 2:])
ax[0].legend(linesL, [r'L', 'S-', 'S0', 'S+', 'X-', 'X0'])
linesX = ax[1].plot(UXlist, rhos_xi[:, 2:])
ax[1].legend(linesX, [r'L', 'S-', 'S0', 'S+', 'X-', 'X0'])

print masses_lambda, masses_xi

plt.show()




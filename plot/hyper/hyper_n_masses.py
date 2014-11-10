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

C = eos.KVOR_mod2()
C.SetHyperConstants(2)
wr = Wrapper(C)
n = np.linspace(0.0, 4.0, 2000)

uX = []
uY = []


C.Csp = 1

C.Cs = 179.56233875157545
C.Co =  87.59963973682763
C.Cr = 100.63642242792484
C.b = 0.007734608051455927
C.c = 0.0003446178665624873
  
C.Cs = 179.56233875157545
C.Co =  141.14599595
C.Cr = 92.11
C.b = 0.004860 
C.c = -0.008312 
  
C.phi_gamma = 3.0
C.phi_z = 3.5

C.omega_c = -1000.0
   
C.Csp = 380.0  
# 
# C.alpha = 0.85
# C.z = 0.65
#         
# C.omega_a = 4.0
# C.omega_f = 0.55
#          
# C.rho_a = 0*150.0
# C.rho_f = 0.35
#      
# C.beta = 0.5
# C.gamma = 8.00
#          
# C.phi_a = -0.5
# C.phi_f = 0.29
#   
# C.d = -2.6

C.alpha = 0.85
C.z = 0.65
C.omega_a = 6.45
C.omega_f = 0.53
          
# C.rho_a = -1500.
# C.rho_f = 0.6
      
C.beta = 0.8
C.gamma = 7.5
          
C.phi_a = -0.85
C.phi_f = 0.28

C.d = -5.5

C.rho_f = 0.75
C.rho_a = 1000
# C.omega_f = 0.6
f0 = 0.27
C.rho_kind = 1
C.rho_power = 2.0
C.omega_c = -15000



suffix = 'Mod/MOD_rhoa_%.2f'%C.rho_a

f0 = 0.27

npoints = 200

print C.f0

C.f0 = f0
wr.solve(f0 = C.f0, iter = 3000)

print C.f0
nmax = 5.0
C.SetHyperConstants(2)
C.set_xs(np.array([1.0, 1., -28., 30. , 30., 30., -15., -15.]))
C.sprime = 0
C.phi_meson = 1
hyper = 1
flist = [f0]#[0.195, 0.26, 0.27, 0.28]
lines = []
labels = []
n_star = []
mlist = []
Vs = []

sc = 1+C.sprime

for f in flist:
    global n_star
#     for _f in np.linspace(C.f0, f):
    wr.solve(f0=f, iter=3000)
#     exit()
    print C.Cs, C.Co, C.Cr, C.b, C.c    
    print 'J = ', eos.K(wr.n0, C)
    print 'K = ',eos.J(wr.n0, C)
    print 'K\' = ',-3*wr.n0*(derivative(lambda z: eos.K(z, C), wr.n0, dx=1e-3, order=3) - 
                   2*eos.K(wr.n0,C)/wr.n0)
    print 'L = ', 3*wr.n0*derivative(lambda z: eos.J(z, C), wr.n0, dx=1e-3)

    lines.append(wr.Psymm(n))
    labels.append(str(C.phi_n(0, f)))
    wr.reset(hyper=hyper, nmax=nmax, npoints=npoints, iter=10)
    
    wr.setDriver()
    n_star, M, R = wr.stars(nmax=nmax, npoints=200)
    mlist.append(M)

rho = []
for i in wr.rho:
    print i
    rho.append(i[sc:]/np.sum(i[sc:]))

rho = np.array(rho)
print rho
plt.plot(wr.n/wr.n0, rho)
plt.show()
print sc
with open(suffix+'hyper_n_%.2f_%i_%i.dat' % (C.f0, C.phi_meson, C.sprime), 'w') as f:
    for i, r in enumerate(rho):
        f.write('%f   ' % (wr.n[i]/wr.n0))
        for _r in r:
            f.write('%f  ' % _r)
        f.write('\n')  

mlist = np.array(mlist)

for j,m in enumerate(mlist):
    with open(suffix+'mass_scaling_%.2f_%i_%i.dat'%(flist[j], C.phi_meson, C.sprime),'w') as f:
        for i, _n in enumerate(n_star):
            print m[i]
            print max(m)
            if i > np.argmax(m):
                break
            f.write('%f  %f\n'%(_n/wr.n0, m[i]))
    MofN = interpolate.interp1d(n_star, m)
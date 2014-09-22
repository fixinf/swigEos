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

C.sprime = 0
C.Csp = 380.0  
  
# C.phi_gamma = 3.0
# C.phi_z = 3.5
# 
# C.omega_c = -1000.0
    
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
#            
C.rho_a = 0*-1000.0
C.rho_f = 0.4
#        
C.alpha=0.85

C.beta = 0.5
C.gamma = 8.00
            
C.phi_a = -0.85
C.phi_f = 0.28
  
C.phi_gamma = 3.0
C.phi_z = 2.0
  
C.d = -5.5
C.omega_c = 0*-1000

C.f0 = 0.27

C.SetHyperConstants(2)
print [i for i in C.X_s]
C.set_xs(np.array([0.0, 0.0, -30.0, 30.0, 30., 30., -18., -18.]))
print [i for i in C.X_s]
wr.solve(f0=C.f0, iter=3000)

wr.reset(hyper=0, npoints=2000)
C.SetHyperConstants(2)

print C.f0
print C.eta_r(C.f0), C.eta_o(C.f0), C.eta_p(C.f0)

fig, ax = plt.subplots(1,2)

rho = []
for r in wr.rho:
    rho.append(r[1:]/np.sum(r[1:]))

rho = np.array(rho)
E, P, N = wr.EPN()
ax[1].plot(E, P)
ax[0].plot(wr.n/wr.n0, rho, wr.n/wr.n0, [0.14 for i in wr.n], wr.n/wr.n0, wr.rho[:, 0])
plt.show()

Jlist = []
f_eq = [0.0]
for i in wr.n:
    f_eq = eos.f_eq(np.array([i/2, i/2]), np.array([f_eq[0]]), 1, C)
    Jlist.append(eos.J(i, C, f_eq[0]))
    
plt.plot(wr.n/wr.n0, Jlist)
plt.plot(wr.n/wr.n0, map(lambda z: C.eta_r(z), wr.rho[:,0]))
# plt.ylim([0.0, 100.0])
plt.show()

wr.setDriver()

N, M, R = wr.stars(npoints=100)
print 'Mmax = ', np.max(M)

frange = np.linspace(0.0, 1.0, 1000)
eta_o = []
eta_r = []
eta_s = []
phi_n = []
eta_p = []
fi = open('scalings_f%.2f.dat'%C.f0, 'w')
for f in frange:
    eta_o.append(C.eta_o(f))
    eta_r.append(C.eta_r(f))
    eta_p.append(C.eta_p(f))
    eta_s.append(C.eta_s(f))
    phi_n.append(C.phi_n(0, f))
    fi.write('%f %f %f %f %f %f \n'%(f, C.eta_s(f), C.eta_o(f), C.eta_r(f), C.eta_p(f), C.phi_n(0,f)))
    
fi.close()
etas = np.array([eta_s, eta_o, eta_r, eta_p, phi_n]).transpose()
lines = plt.plot(frange, etas)
plt.legend(lines,['s', 'o', 'r', 'p', 'n'], loc=0)
plt.show()

eta_o = []
eta_r = []
eta_s = []
phi_n = []
eta_p = []
print C.f0
pause(5)
fi = open('scalings_NSf%.2f.dat'%C.f0, 'w')
for i, f in enumerate(wr.rho[:,0]):
    eta_o.append(C.eta_o(f))
    eta_r.append(C.eta_r(f))
    eta_p.append(C.eta_p(f))
    eta_s.append(C.eta_s(f))
    phi_n.append(C.phi_n(0, f))
    fi.write('%f %f %f %f %f %f \n'%(wr.n[i]/wr.n0, C.eta_s(f), C.eta_o(f), C.eta_r(f), C.eta_p(f), C.phi_n(0,f)))
    
fi.close()
etas = np.array([eta_s, eta_o, eta_r, eta_p, phi_n]).transpose()
lines = plt.plot(wr.n/wr.n0, etas)
plt.legend(lines,['s', 'o', 'r', 'p', 'n'], loc=0)
plt.show()

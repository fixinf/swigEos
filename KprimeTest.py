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

C.Csp = 1

C.Cs = 179.56233875157545
C.Co =  141.14599595
C.Cr = 92.11
C.b = 0.004860 
C.c = -0.008312 
  
# C.alpha = 0.85
# C.z = 0.65
# 
# C.omega_a = 6.45
# C.omega_f = 0.53
#         
# C.rho_a = 0*150.0
# C.rho_f = 0.35
#     
# C.beta = 0.5
# C.gamma = 8.00
#---------------------     
C.phi_a = -0.85
C.phi_f = 0.28
  
C.d = -5.5

C.alpha = 0.85
C.z = 0.65

C.phi_gamma = 3.0
C.phi_z = 3.5
     
C.sprime = 0
C.Csp = 380.0
  
C.f0 = 0.26
 
C.beta = 1.0
C.gamma = 14.0

f0 = 0.26
  
wr.solve(f0=C.f0)
C.SetHyperConstants(2)
npoints = 1000
dn = 1e-3

print -3*wr.n0*((eos.K(wr.n0+dn, C) - eos.K(wr.n0-dn, C))/(2*dn) - 2*eos.K(wr.n0, C)/wr.n0) 

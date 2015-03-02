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

cAPR = 1.0/0.16
nAPR=cAPR*np.array([0.0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.32, 0.4, 0.48, 0.56, 0.64, 0.8, 0.96])
APR = np.array([0.0, -6.48, -12.13, -15.04, -16.00, -15.09, -12.88, -5.03, 2.13, 15.46, 
             34.39, 58.35, 121.25, 204.02])

nrange = np.linspace(nAPR[0], nAPR[-1], 50)
iAPR = interpolate.interp1d(nAPR, APR, kind='cubic')

with open('APR.dat', 'w') as f:
    for i, n in enumerate(nrange):
        f.write('%f %f \n'% (n, iAPR(n)))
# exit()    

APR_N = np.array([0., 4.45, 6.45 , 9.65, 13.29, 17.94, 22.92, 27.49, 38.82, 54.95, 75.13, 99.75, 127.58, 205.34, 305.87])
nAPR_N = cAPR* np.array([0., 0.02, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24,  0.32,  0.4, 0.48, 0.56, 0.64, 0.8, 0.96]) 

nrange = np.linspace(nAPR_N[0], nAPR_N[-1], 50)
iAPR_N = interpolate.interp1d(nAPR_N, APR_N, kind='cubic')

with open('APR_N.dat', 'w') as f:
    for i, n in enumerate(nAPR_N):
        f.write('%f %f \n'% (n, iAPR_N(n)))
exit()


uX = []
uY = []
C.Csp = 1

C.Cs = 227.825938
C.Co =  134.882548
C.Cr = 93.199060
C.b = 0.005592
C.c = 0.012123

C.alpha = 0.85
C.z = 0.65
C.omega_a = 6.45
C.omega_f = 0.53
          
C.beta = 0.8
C.gamma = 7.5
          
C.phi_a = -0.85
C.phi_f = 0.28

C.d = -5.5

C.rho_f = 0.75
C.rho_a = 1000
f0 = 0.27
C.rho_kind = 1
C.rho_power = 2.0
C.omega_c = -15000

suffix = 'Mod/KVOR_rhoa_%.2f'%C.rho_a

C.SetHyperConstants(2)

npoints = 1000

with open('../klahnUpper2', 'r') as f:
    for line in f:
        x,y = line.split()
        uX.append(float(x)/0.16)
        uY.append(float(y))


fig, ax = plt.subplots()

print C.f0
C.f0 = f0
wr.solve(f0 = C.f0, iter = 3000)
print C.f0

C.SetHyperConstants(2)

flist = [f0]

ax.plot(nAPR, APR)
ax.set_ylim([-30, 300])
Es = wr.Esymm(n)
Es /= n
Es -= C.M[0]
Es *= wr.m_pi
print Es
ax.plot(n/wr.n0, Es)
plt.show()
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
n = np.linspace(0.0, 4.0, 5000)

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

C.Cs = 179.56233875157545
C.Co =  87.59963973682763
C.Cr = 100.63642242792484
C.b = 0.007734608051455927
C.c = 0.0003446178665624873

C.sprime = 0
C.Csp = 380.0  
# 

wr.solve()
fig, ax = plt.subplots(1, 2)
plt.plot(n[1:]/wr.n0, wr.Psymm(n))
# C.phi_f = 0.35
# C.phi_a = -8.0
C.d = -1.5
wr.solve()
plt.plot(n[1:]/wr.n0, wr.Psymm(n))
plt.ylim([-100, 1000])
plt.show()
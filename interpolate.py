#!/usr/bin/python

import eosWrap as eos
from Wrapper import Wrapper
from scipy import interpolate
from pylab import plot, show
import numpy as np

C = eos.KVOR()
C.Csp = 1.0
wr = Wrapper(C)
wr.reset(nmax=5.5, npoints = 1000)
E, P, N = wr.EPN()
N = N/wr.n0
e = []
p = []
n = []

with open("crust.dat", 'r') as f:
    for line in f:
        _e, _p, _n = line.split()
        e.append(float(_e))
        p.append(float(_p))
        n.append(float(_n))

print P
# 
# plist = np.append(p[:],P[:])
# elist = np.append(e[:],E[:])
# nlist = np.append(n[:],N[:])

# 
plist = P[:]
elist = E[:]
nlist = N[:]

iP = interpolate.interp1d(nlist, plist)
iE = interpolate.interp1d(nlist, elist)

print nlist

nmax = 11.0
print nlist[0]
print iP(0.5)
iN = np.linspace(0.0, nmax, 100)


# plot(iN, map(iE, iN))
# show()

dr = eos.KVDriver()
dr.set(E, P, N*wr.n0)
nstar = np.linspace(0.5, 5.0, 20)
marr = np.array(map(lambda z: eos.star_crust(z, 2, dr), nstar))
dr.set(iE(iN), iP(iN), iN*wr.n0)
marr2 = np.array(map(lambda z: eos.star_crust(z, 2, dr), nstar))
print marr[:,0]
print marr2[:,0]

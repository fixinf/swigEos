#!/usr/bin/python
import eosWrap as eos
from Wrapper import Wrapper
from scipy import interpolate
from pylab import plot, show
import numpy as np

n_cut_crust = 0.7
n_cut_eos = 0.8

C = eos.KVOR()
wr = Wrapper(C)
wr.reset(nmin = n_cut_eos*wr.n0, npoints = 2000)


E, P, N = wr.EPN()
N = N/wr.n0
e = []
p = []
n = []



with open("crust.dat", 'r') as f:
	for line in f:
		# print line
		_e, _p, _n = line.split()
		if float(_n) < n_cut_crust:
			e.append(float(_e))
			p.append(float(_p))
			n.append(float(_n))
n_eos = 3
plist = np.append(p[:],P[:n_eos])
elist = np.append(e[:],E[:n_eos])
nlist = np.append(n[:],N[:n_eos])

iP = interpolate.interp1d(nlist, plist, kind='quadratic')
iE = interpolate.interp1d(nlist, elist, kind='quadratic')

nmax = 7.0
print nlist[0]
print iP(0.5)
gamma = 1./4.
iN = np.linspace(1e-10**gamma, n_cut_eos**gamma, 2000)
iN = iN**(1./gamma)

crust_p = np.array(map(iP, iN))
crust_e = np.array(map(iE, iN))

finalE = np.append(crust_e, E)
finalP = np.append(crust_p, P)
finalN = np.append(iN, N)

# plot(finalN, finalP, n, p, N, P)
# show()

dr = eos.KVDriver()
dr.set(finalE, finalP, finalN*wr.n0)
nstar = np.linspace(2e-1, 4.0, 100)
res = eos.star_crust(1.2905, 3, dr)
print res
print res[2]*wr.m_pi**4 * C.M[0]
exit()
k = 0
def star_print(z):
	global k
	print k
	k+=1
	return eos.star_crust(z, 2, dr)

MR = np.array(map(star_print, nstar))
plot(MR[:,1], MR[:,0])
print max(MR[:, 0])
show()

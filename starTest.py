import numpy as np
import eosWrap as eos
import Models
from Wrapper import Wrapper
from matplotlib import pyplot as plt
from pylab import sqrt, pi
C = Models.KVOR()
wr = Wrapper(C)
# wr.reset(nmin = 0., npoints=1000)
# 
# dr = eos.KVDriver()
# dr.set(wr.E, wr.P, wr.n)
# 
# dr2 = eos.KVDriver()
# E,P,N = wr.EPN()
# dr2.set(E, P, N)


nstar, M, R, Mg1, Mg2 = wr.stars_crust()
plt.plot(Mg1, M, Mg2, M)
plt.show()
exit()
st = eos.star_crust2(1.39, 3, dr, 0.8*wr.n0)
print st
print st[2] * (931.5/135.)
print dr.nSize
print dr.lastNstar
lastN = dr.getLastN(dr.nSize)[:-1]
lastR = dr.getLastR(dr.nSize)[:-1]
lastM = dr.getLastM(dr.nSize)[:-1]
print lastN
print lastR
print lastM
# plt.plot(lastR, lastN, lastR, lastM)
# plt.show()
dx = lastR[1] - lastR[0]
print 'dx=' , dx
print np.trapz(lastN, dx=dx)
grav_mult = []

for i, r in enumerate(lastR):
    grav_mult.append( r**2 / sqrt(1 - 2 * 1.4677 * lastM[i] / r))

grav_mult = np.array(grav_mult)
print grav_mult
res = np.multiply(lastN, grav_mult)

print (0.0004898007281478712)*np.trapz(res, dx=dx)

# print eos.star2(3., 2, dr)
# st2 = eos.star_crust(1.39, 3, dr2, 0.8*wr.n0)
# print st2
# st = eos.star_crust2(1.39, 3, dr, 0.8*wr.n0)
# print st
# print st[2] * (931.5/135.)
# print st2[2] * wr.m_pi**3*931.5
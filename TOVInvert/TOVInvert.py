import matplotlib
from scipy.interpolate.interpolate import interp1d
matplotlib.use('Qt4Agg')
from matplotlib import pyplot as plt
import numpy as np
from Wrapper import Wrapper
import Models
from tabulate import tabulate
from scipy import integrate
from scipy import interpolate
from pylab import pi, sqrt
from Invert import EnthalpySolver

wr = Wrapper(Models.KVOR())
# n, m, r, mb1, mb2 = wr.stars_crust()
# tab = np.array([m, r]).transpose()
# table = tabulate(tab, [], tablefmt='plain') 
# with open('mr.dat', 'w') as f:
#     f.write(table)

n_cut = 0.8*wr.n0
nmax = 5.38721*wr.n0
nrange = np.linspace(n_cut, nmax, 200)
eps_n = wr.Eneutr(nrange)
porig_n = wr.P_N(nrange)/wr.const/wr.m_pi**4

wr.reset(nmin=n_cut, nmax=nmax, npoints=200)
eps = wr.E
porig = wr.P

s = EnthalpySolver(E=eps, P=porig, nhpoints=100)
_m, _r, _p, _e = s.integrateOut(porig[-1])
mmax = _m[-1]
rmax = _r[-1]
print '###########################################'
# print s.integrateIn(mmax, rmax, s.iHofP(porig[-2]))

# m_new = 1.92946
# r_new = 10.596
 
# target_e = 23.7753   
# target_p = 7.56932  


m_new = mmax
r_new = rmax



print 'M = ', mmax, 'R = ', rmax
_mcore, _rcore = s.integrateIn(m_new, r_new, s.iHofP(porig[-1]))
mcore = _mcore[-1]
rcore = _rcore[-1]

lines = plt.plot(_rcore, _mcore, _r, _m)
plt.legend(lines, ['in', 'out'])
plt.show()

mcore2 = _m[1]
rcore2 = _r[1]

print 'M_core = ', mcore, 'R_core =' , rcore
print 'E_int = ' , mcore/(rcore**3 / 3 * s.E_const)
print 'M_core2 = ' , mcore2, 'R_core2 =' , rcore2
print 'E_int2 = ' , mcore2/(rcore2**3 / 3 * s.E_const)
print 'P = ', porig[-1], 'E = ', eps[-1]
print s.predictPE(mcore, rcore, eps[-1], porig[-1])

exit()


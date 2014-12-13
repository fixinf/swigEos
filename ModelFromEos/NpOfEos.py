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
import eosWrap as eos
from scipy.misc import derivative

C = Models.waleckaMatsui()
wr = Wrapper(C)
wr.reset(npoints=100)


E = wr.E[3]
P = wr.P[3]
n = wr.n[3]

cAPR = 1.0/0.16
nAPR=cAPR*np.array([-1e-2, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.32, 0.4, 0.48, 0.56, 0.64, 0.8, 0.96])
APR = np.array([0.0, -6.48, -12.13, -15.04, -16.00, -15.09, -12.88, -5.03, 2.13, 15.46, 
             34.39, 58.35, 121.25, 204.02])

APR_N = np.array([0., 4.45, 6.45 , 9.65, 13.29, 17.94, 22.92, 27.49, 38.82, 54.95, 75.13, 99.75, 127.58, 205.34, 305.87])
nAPR_N = cAPR* np.array([0., 0.02, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24,  0.32,  0.4, 0.48, 0.56, 0.64, 0.8, 0.96]) 

aprEps = nAPR*wr.n0*APR/wr.m_pi + C.M[0]*nAPR*wr.n0
i_apr_eps = interp1d(nAPR*wr.n0, aprEps, kind='cubic')
print nAPR*wr.n0
i_apr_P = lambda z: z*derivative(i_apr_eps, z, dx=1e-3) - i_apr_eps(z)
# nrange = np.linspace(1e-2, (0.95/0.16)*wr.n0, 100, endpoint=False)
# tab = np.array([nrange/wr.n0, map(i_apr_eps, nrange), map(i_apr_P, nrange)]).transpose()
# table = tabulate(tab, ['n/n_0', 'E [m_pi^4]', 'P [m_pi^4]'], tablefmt='plain')
# 
# with open('APR_N.dat', 'w') as f:
#     f.write(table)

print  i_apr_eps(wr.n0), i_apr_P(wr.n0)

# E = float(i_apr_eps(wr.n0))
# P = float(i_apr_P(wr.n0))
# n = wr.n0

C.Co = C.Co
C.Cr = C.Cr
print 'C_o^2 = ', C.Co
print eos.solveF(n, E, P, np.array([0.1, 0.1,]), 2, C)
print n, E, P
print wr.rho[3]


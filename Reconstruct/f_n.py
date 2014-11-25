import matplotlib
from tabulate import tabulate
from gtk.keysyms import PesetaSign
matplotlib.use('QT4Agg')
from math import pi
from numpy import sqrt, arange, array
from pylab import *
from scipy import optimize as opt
from scipy.interpolate import interp1d
from scipy.integrate.quadpack import quad
from scipy.misc.common import derivative

import eosWrap as eos
import Models
from Wrapper import Wrapper
from scipy.integrate.odepack import odeint

m_pi=135.0
m_N = 939.0/m_pi
m_om = 782.0/m_pi
Co = 54.6041
Cs = 164.462
Cr = 121.69
b = 0.0202832
c = 0.0471633
d = 0.0

C = eos.Walecka()
C.Cs = Cs
C.Co = Co
C.Cr = Cr
C.b = b
C.c = c

wr = Wrapper(C)

def pf(n):
    return (3 * pi**2 * n)**(1.0/3.0)

def K(n, m):
    p_f = pf(n)
    res = p_f*sqrt(m*m + p_f*p_f)*(m*m + 2*p_f*p_f);
    if (m > 0.0):
        res -= m**4*arcsinh(p_f/m);
    res /= (8*pi*pi);
    return res

def U1(f):
    return m_N**4*(b*f**3/3 + c * f**4 / 4) 

def E(n, f, eta_o=None, U=None):
    res = m_N**4 * f**2 / (2 * Cs)
    om = Co * n**2/(2*m_N**2)
    if eta_o is not None:
        om /= eta_o(f)
#         print '!'
    res += om
    res += 2*K(n/2, m_N*(1-f))
    if U is not None:
        res += U(f)
    else:
        res += U1(f)
    return res*135

def EN(n, f, eta_o=None, U=None):
    res = m_N**4 * f**2 / (2 * Cs)
    om = Co * n**2/(2*m_N**2)
    if eta_o is not None:
        om /= eta_o(f)
#         print '!'
    res += om
    res += Cr*n**2/(8*m_N**2)
    res += K(n, m_N*(1-f))
    if U is not None:
        res += U(f)
    else:
        res += U1(f)
    return res*135


def E2(n, f, o, U = None):
    res = E(n, f, U=U)
    res -= Co * n**2/(2*m_N**2)
    g_o = sqrt(Co)*m_om/m_N
    res += g_o * o * n
    res -= m_om**2 * o**2 / 2
    res += d*o**4
    return res

def omega(n, f):
    dx = 1e-5
    fun = lambda z: (E2(n, f, z + dx) - E2(n, f, z))/dx
    return opt.newton(fun, 0.5, maxiter=200)

def fun_f_eq2(f, n):
    dx = 1e-6
    o = omega(n,f)
    return (E2(n, f+dx, o) - E2(n, f, o))/dx

def f_eq2(n):
    return opt.newton(lambda z: fun_f_eq2(z, n), 0.1, maxiter = 200)

def f_eq(n):
    dx = 1e-6
    fun = lambda z: (E(n, z+dx) - E(n, z))/dx
    return opt.newton(fun, 0.15, maxiter=200)

def f_eq_N(n):
    dx = 1e-6
    fun = lambda z: (EN(n, z+dx) - EN(n, z))/dx
    return opt.newton(fun, 0.15, maxiter=200)
# 
# print (E(0.5, 0.15)/0.5 - m_N)*135.0
# print omega(0.5, 0.15)
# print ((E2(0.5, 0.15, omega(0.5, 0.15))/0.5 - m_N)*135.0)
# print f_eq(1.0)
nrange = linspace(0.0, 4.0, 100)
n_inter = linspace(0.01, 4.5, 100)

flist = map(lambda z: f_eq(z), nrange)

pflistN = map(lambda z: pf(z), nrange)
pflistS = map(lambda z: pf(z), nrange/2)

my_fpS = interp1d(pflistS, flist, kind='cubic')
my_fpN = interp1d(pflistN, flist, kind='cubic')

print my_fpS(3.5)

def func_integr(z, f):
#     print 'z=', z
    return (z*z/pi**2)*sqrt(z**2 + m_N**2 * (1 - f(z))**2)



def _E_N(n):
    res = Co*n**2/(2*m_N**2)
    int = quad(lambda z: func_integr(z, my_fpN), 0.0, pf(n))[0]
    res += Cr*n**2/(8*m_N**2)
#     print 'int = ', int
    res += int
    return res*135



########################Constructing f(n)#############################
# f0 = 0.15
# n0 = 0.5
# K0 = 235.0
# E0 = -16.0
# sq = sqrt(pf(n0/2)**2 + m_N**2*(1-f0)**2)
# C_om = m_N**2*((E0/135.0 + m_N)/n0 - (1.0/n0) * sq)
# print C_om
# C1 = (sq*(K0 - 9*n0**2 * C_om/(m_N**2)) - (3*pi**2)**(2.0/3.0)*(n0/2)**(1.0/3.0))/(-m_N**2 * (1-f0))
# C2 = 0.25*(E0/135.0 - m_N) * n0 + 0.5*n0*sq
# print C1, C2
# def  func_optimize(alpha):
#     beta = -(2/n0**2)*(-6*f0 + 2 * alpha * n0 + C1*n0)
#     gamma = (4/n0**3)*(-4 * f0 + alpha*n0 + C1*n0)
#     print 'beta=', beta, 'gamma=', gamma
#     my_f = lambda p: alpha*p + beta*p**2 + gamma*p**3
#     res = quad(lambda z: z*z*sqrt(z*z + m_N**2*(1-my_f(z))), 0, pf(n0/2))[0]
#     return res - C2
# 
# print func_optimize(0.0)
# show()

#######################Fitting to Landau params######################
C2 = Models.waleckaMatsui()
# C2 =  Models.KVOR()
wr = Wrapper(C2)
# 
# lomb_n = []
# lomb_f0 = []
# with open('lomb_params.dat','r') as f:
#     for line in f:
#         k, f0 = line.split()
#         print k, f0, float(k)
#         lomb_n.append( (float(k)*197.33/135.)**3 * 2 /(3*pi**2))
#         lomb_f0.append(float(f0))
# lomb_n = np.array(lomb_n)
# print lomb_n
# 
# 
# nrange2 = lomb_n#np.linspace(lomb_n[0], lomb_n[-1], 100)
# nrange = lomb_n[1:-1]#np.linspace(lomb_n[0]-0.01, lomb_n[-1]+0.5, 100)
# F0 = wr.f0(nrange2, multiply=False)

# lomb_f0 = np.array(lomb_f0)
# print F0
# plt.plot(lomb_n/wr.n0, lomb_f0, nrange2/wr.n0, F0)
# plt.show()

nrange2 = np.linspace(0., 4.5, 100)
nrange = np.linspace(0., 2., 100)
F0 = wr.f0(nrange2, multiply=False)


iF0 = interp1d(nrange2, F0)
# iF0 = interp1d(lomb_n, lomb_f0)

Co = C2.Co*0.7

def func_f_solve(f, n):
    f = float(f)
    print f, n
    pf = eos.p_f(n/2)
    mn = C.M[0]
    meff = mn * C.phi_n(0, f)
    res = -sqrt(pf**2 + meff**2)/(mn**2 * C.phi_n(0, f))
    res *= (iF0(n) - Co/mn**2)
    return res

# print func_f_solve(0.2, 0.2)

flist = odeint(func_f_solve, 0., nrange)

f2 = []
f= 0.
for n in nrange:
    f, =eos.f_eq(np.array([n/2, n/2]), np.array([f]), 1, C2)
    f2.append(f)

dn = nrange[1] - nrange[0]
df2 = np.diff(f2) / dn

func = []
ebind = []
for i, n in enumerate(nrange):
    func.append(func_f_solve(f2[i], n))
    ebind.append(eos.EBind(np.array([f2[i], n/2, n/2]), C2))

pflistS = np.array(map(lambda z: pf(z), nrange/2))
print pflistS.shape, flist.shape
print pflistS[0], pflistS[-1]
#my_fpS = interp1d(pflistS, flist[:,0], kind='cubic')
my_fpS = interp1d(pflistS, flist[:,0], kind='cubic', bounds_error=False, fill_value=pflistS[0])
# my_fpS = interp1d(pflistS, f2, kind='cubic', bounds_error=False, fill_value=pflistS[0])

def _E(n, fp):
    res = Co*n**2/(2*C.M[0]**2)
#     print pf(n/2)
    int = 2*quad(lambda z:func_integr(z, fp), 0.0, pf(n/2))[0]
#     print 'int = ', int
    res += int
    return res

##########################################################################3
### TEST E(n), P(n) -> f(n)
cAPR = 1.0/0.16
nAPR=cAPR*np.array([-1e-2, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.32, 0.4, 0.48, 0.56, 0.64, 0.8, 0.96])
APR = np.array([0.0, -6.48, -12.13, -15.04, -16.00, -15.09, -12.88, -5.03, 2.13, 15.46, 
             34.39, 58.35, 121.25, 204.02])

aprEps = nAPR*wr.n0*APR/wr.m_pi + C.M[0]*nAPR*wr.n0
i_apr_eps = interp1d(nAPR*wr.n0, aprEps, kind='cubic')
print nAPR*wr.n0
i_apr_P = lambda z: z*derivative(i_apr_eps, z, dx=1e-3) - i_apr_eps(z)

print  i_apr_eps(wr.n0), i_apr_P(wr.n0)

wrEs, wrF= wr.Esymm(nrange, ret_f=1)
pEs = wr.Psymm(nrange) / wr.const / wr.m_pi**4



# plt.plot(i_apr_eps(nrange), i_apr_P(nrange), wrEs, pEs)
# plt.show()

wrEs = i_apr_eps(nrange)
pEs = i_apr_P(nrange)

fEs = 1 - sqrt((wrEs + pEs - Co * nrange**2 / C.M[0]**2)**2 / (nrange**2) - pflistS**2)/C.M[0]
fEs = np.nan_to_num(fEs)
fEs[0] = 0.
print fEs
print wrF
lines = plot(nrange, fEs, nrange, wrF)
legend(lines, ['my', 'orig'])
show()
ifEs = interp1d(pflistS, fEs, kind='cubic', bounds_error=False, fill_value=pflistS[0])
print fEs
eInt = map(lambda z: _E(z, ifEs), nrange)
print wrEs
print eInt
lines = plot(nrange, eInt, nrange, wrEs)
legend(lines, ['my', 'orig'])
show()



ebind_new = map(lambda z: _E(z, my_fpS)/z - m_N*135, nrange)

fig, ax = plt.subplots(1, 3)
lines_f = ax[0].plot(nrange, flist, nrange, f2)
# ax[0].plot(nrange, f2)

ax[0].legend(lines_f, ['New f', 'Orig'], loc=0)

lines_func = ax[1].plot(nrange, func, nrange[1:], df2)
# ax[1].plot()
ax[1].legend(lines_func, ['df/dn new', 'df/dn orig'])

lines_e = ax[2].plot(nrange, ebind, nrange, ebind_new)
# ax[2].plot(nrange2, ebind_new)
ax[2].legend(lines_e, ['Ebind new', 'Ebind orig'])
plt.show() 

tab = np.array([nrange, ebind, ebind_new, f2, flist]).transpose()

table = tabulate(tab, ['n/n_0','ebind_old', 'ebind', 'f_old', 'f'], tablefmt='plain')

with open('f_nCo=%.2f.dat'%(Co), 'w') as f:
    f.write(table)

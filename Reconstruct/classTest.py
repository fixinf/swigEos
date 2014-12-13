import matplotlib
from scipy import interpolate
from scipy.misc.common import derivative
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import numpy as np
from Wrapper import Wrapper, FromEos, EosConstructor
import Models
C = Models.waleckaMatsui()
wr = Wrapper(C)
n = np.linspace(0., 8*wr.n0, 80, endpoint=0)
from scipy.interpolate import interp1d
import eosWrap as eos
from pylab import sqrt
#################### TEST ES ######################
E, fs0 = wr.Esymm(n, ret_f=1)
P = wr.Psymm(n)/wr.const/wr.m_pi**4
f_eq = 0.
for i,_n in enumerate(n):
    f_eq, = eos.f_eq(np.array([_n/2, _n/2]), np.array([f_eq]), 1, C) 
    print E[i] + P[i], _n * eos.mu(np.array([f_eq, _n/2, _n/2]), 1, C), _n*(
          C.Co * _n / C.M[0]**2 + sqrt(eos.p_f(_n/2)**2 + C.M[0]**2*(1-f_eq)**2))
    f_eq_new = (E[i] + P[i] - C.Co * _n**2 / C.M[0]**2)


Eos1 = FromEos()
print Eos1.mn, C.M[0]
# exit()
# Eos1.setScaling(lambda z: (1 + 1.65*0.2)/(1 + 1.65*z),
#                 lambda z: 1.)
Eos1.setConstants(Co=C.Co)

Eos2 = EosConstructor()
Eos2.setConstants(Co=Eos1.Co)
fs = Eos1.FsFromEos(E, P, n)

plt.plot(n/wr.n0, fs, n/wr.n0, fs0)
plt.show()
exit()

Eos2.setFs(fs0, n)
Es2 = map(lambda z: Eos1.ES(z), n)

lines = plt.plot(n, E, n, Es2)
plt.legend(lines, ['old', 'new'], loc=0)
plt.show()

#################### TEST EN ######################

E, fn0 = wr.Eneutr(n, ret_f=1)
P = wr.P_N(n)/wr.const/wr.m_pi**4

Eos1 = FromEos()
# Eos1.setScaling(lambda z: (1 + 1.65*0.2)/(1 + 1.65*z),
#                 lambda z: 1.)
Eos1.setConstants(Co=C.Co)

Eos2 = EosConstructor()

fn = Eos1.FnFromEos(E, P, n)
# plt.plot(n/wr.n0, fn, n/wr.n0, fn0)
# plt.show()
# exit()
Eos2.setFn(fn0, n)
Eos2.setConstants(Co=C.Co, Cr=C.Cr)
# En2 = map(lambda z: Eos2.EN(z), n)
# lines = plt.plot(n, E, n, En2)
# plt.legend(lines, ['old', 'new'], loc=0)
# plt.show()

#################### TEST F(N,X) EXPANSION ########
# iFs0 = interpolate.interp1d(n, fs0, kind='cubic')
# iFn0 = interpolate.interp1d(n, fn0, kind='cubic')
# iFx0 = lambda z, x: iFs0(z) + (iFn0(z) - iFs0(z)) * (1 - 2*x)**2
# wr.reset(nmax=3.5)
# res = []
# for _i, _n in enumerate(wr.n):
#     res.append(iFx0(_n, wr.rho[_i, 2]/_n))
#     
# print res
# lines=plt.plot(wr.n, res, wr.n, wr.rho[:,0], wr.n, iFs0(wr.n))
# plt.legend(lines, ['inter','old','symm'])
# plt.show()

#################### TEST APR RECONSTRUCTION###########
cAPR = 1.0/0.16
nAPR=cAPR*np.array([-1e-2, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.32, 0.4, 0.48, 0.56, 0.64, 0.8, 0.96])
APR = np.array([0.0, -6.48, -12.13, -15.04, -16.00, -15.09, -12.88, -5.03, 2.13, 15.46, 
             34.39, 58.35, 121.25, 204.02])

APR_N = np.array([0., 4.45, 6.45 , 9.65, 13.29, 17.94, 22.92, 27.49, 38.82, 54.95, 75.13, 99.75, 127.58, 205.34, 305.87])
nAPR_N = cAPR* np.array([-1e-2, 0.02, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24,  0.32,  0.4, 0.48, 0.56, 0.64, 0.8, 0.96]) 

aprEps = nAPR*wr.n0*APR/wr.m_pi + C.M[0]*nAPR*wr.n0
aprEpsN = nAPR_N*wr.n0*APR_N/wr.m_pi + C.M[0]*nAPR_N*wr.n0

i_apr_eps = interp1d(nAPR*wr.n0, aprEps, kind='cubic')
i_apr_P = lambda z: z*derivative(i_apr_eps, z, dx=1e-3) - i_apr_eps(z)

i_apr_epsN = interp1d(nAPR_N*wr.n0, aprEpsN, kind='cubic')
i_apr_PN = lambda z: z*derivative(i_apr_epsN, z, dx=1e-3) - i_apr_epsN(z)


nrange = np.linspace(0.0, 5.9*wr.n0, 100)
E = i_apr_eps(nrange)
P = i_apr_P(nrange)

EN = i_apr_epsN(nrange)
PN = i_apr_PN(nrange)

rAPR = FromEos()
Co = 100
Cr = 100
rAPR.setConstants(Co=Co, Cr=Cr)
fsAPR = rAPR.FsFromEos(E, P, nrange)
fnAPR = rAPR.FnFromEos(EN, PN, nrange)

x_apr = [0.]
eApr = [0.]
with open('AkmalBetaEqSl.csv', 'r') as f:
    for _i, line in enumerate(f):
        if _i > 0:
            x, e = line.split(',')
            x_apr.append(float(x))
            eApr.append(float(e))
x_apr = np.array(x_apr)
eApr = np.array(eApr)

ieApr = interpolate.UnivariateSpline(x_apr, x_apr*eApr*wr.n0/0.16/135 + C.M[0]*x_apr*wr.n0/0.16)
ipApr = lambda z: z*derivative(ieApr, z, dx=1e-3) - ieApr(z)
plt.plot(x_apr[2:]/0.16, np.diff(ipApr(x_apr[1:]))/np.diff(ieApr(x_apr[1:])))
plt.show()

nrange = x_apr/0.16*wr.n0
        
ENS = ieApr(nrange/wr.n0*0.16)
PNS = ipApr(nrange/wr.n0*0.16)

fnsApr = rAPR.FnFromEos(ENS, PNS, nrange)
plt.plot(nrange/wr.n0, fnsApr)
ens=map( rAPR.EN, nrange[1:])
plt.show()
plt.plot(nrange[1:]/wr.n0, ens, nrange/wr.n0, ENS)
plt.show()
i_ens = interpolate.interp1d(nrange, ens)
pns = lambda z: z * derivative(lambda x: i_ens(x),z, dx=1e-3) - i_ens(z)
print ens - ENS
print pns
plt.plot(ens[1:], pns(nrange), ENS, PNS)
plt.show()
# apr_vsN = np.diff(P)/np.diff(E)
# plt.plot(EN, PN, EN, EN)
# plt.show()
# plt.plot(nrange, fsAPR, nrange, fnAPR)
# plt.show()
# plt.plot(nrange, map(rAPR.ES, nrange), nrange, E)
# plt.show()

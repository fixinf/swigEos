import eosWrap as eos
import matplotlib
from scipy.misc.common import derivative
from tabulate import tabulate
matplotlib.use('QT4Agg')
from Wrapper import Wrapper
from scipy.integrate.quadpack import quad
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import Models

# C = Models.KVOR()
fname_0='KVOR/f0.dat'
fname_1='KVOR/f1.dat'
C = Models.myMod()
# C = eos.Walecka()
C.SetHyperConstants(2)
C.Csp = 1.



# C.Cs = 266.9
# C.Co = 195.7
# C.b = 0.
# C.c = 0.

wr = Wrapper(C)

def I1(n, m):
    p = eos.p_f(n)
    res = sqrt(m**2 + p**2)*(-3*m**2 * p/8. + p**3 / 4.)
    res += 3 * log(p/m + sqrt(1 + p**2/m**2)) / 8.
    res /= (pi**2)
    return res

def I1num(n, m):
    p = eos.p_f(n)
    res = quad(lambda z: z**4/(z**2 + m**2)**(3./2), 0, p)[0]
    res /= pi**2
    return res


print I1(1., 1.)
print I1num(1., 1.)

def I2(n, m):
    p = eos.p_f(n)
    res = p * sqrt(p**2 + m**2)*(m**2 + 2*p**2)
    res -= m**4 * log(p/m + sqrt(1 + p**2 / m**2))
    res /= 8* pi **2
    return res

def I2num(n, m):
    p = eos.p_f(n)
    res = quad(lambda z: z**2*sqrt(z**2 + m**2), 0, p)[0]
    res /= pi**2
    return res

def I3(n, m):
    p = eos.p_f(n)
    res = p * sqrt(p**2 + m**2)
    res -= m**2 * log(p/m + sqrt(1 + p**2 / m**2))
    res /= 2 * pi **2
    return res

def I3num(n, m):
    p = eos.p_f(n)
    res = quad(lambda z: z**2/sqrt(z**2 + m**2), 0, p)[0]
    res /= pi**2
    return res

def I4num(n, m):
    p = eos.p_f(n)
    res = quad(lambda z: z**2/(z**2 + m**2)**(3./2), 0, p)[0]
    res /= pi**2
    return res

print I3(1., 1.)
print I3num(1., 1.)

print I2(1., 1.)

print I2num(1., 1.)

print I4num(1., 1.)
# exit()
f_f0 = 0.

denomlist = []
flist = []
numlist = []
def f0(n):
    global f_f0
    mn = C.M[0]
    pf = eos.p_f(n/2)
    f_f0, = eos.f_eq(np.array([n/2, n/2]), np.array([f_f0]), 1, C)
    f = f_f0
    flist.append(f)
    meff = mn * C.phi_n(0, f) 
    res = C.Co/(mn**2 * C.eta_o(f))
    dPhi = derivative(lambda z: C.phi_n(0, z), f, dx=1e-3)
    dEtaO = derivative(lambda z: C.eta_o(z), f, dx=1e-3)
#     print 'dEtaO = ', dEtaO
    U = lambda z: C.U(z) + mn**4 * z**2 * (C.eta_s(z))/ (2* C.Cs)
    d2U = derivative(lambda z: U(z), f, dx=1e-3, n=2)
    d2EtaO = derivative(lambda z: C.eta_o(z), f, dx=1e-3, n=2)
    d2Phi = derivative(lambda z: C.phi_n(0, z), f, dx=1e-3, n=2)
#     print n, d2Phi, d2EtaO, 2* C.Cs/mn**2 * d2Phi * C.phi_n(0,f)* I3num(n/2, meff)
    eta_contrib = derivative(lambda z: 1./C.eta_o(z), f, dx=1e-3, n=2)#2*dEtaO**2 / C.eta_o(f)**3 - d2EtaO / C.eta_o(f)**2
#     print 'n= %f, f = %f, eta_contrib = %.12f'%(n, f, eta_contrib)
    res2 = -C.Cs / mn**4 #* C.phi_n(0, f) * dPhi/ (sqrt(pf**2 + meff**2) * mn**2)
    num = (C.Co*n*dEtaO/(mn**2 * C.eta_o(f)**2) - mn**2 * C.phi_n(0, f) * dPhi / sqrt(pf**2 + meff**2))**2
    
    res2 *= num
#     res2 *= - mn**2 * C.phi_n(0, f) * dPhi / sqrt(pf**2 + meff**2)
    denom = (2 * C.Cs / mn**2 * dPhi**2 * I1num(n/2, meff) + d2U * C.Cs/mn**4 + 
        C.Cs * C.Co * n**2 * eta_contrib / (2*mn**6) + 
        2* C.Cs/mn**2 * d2Phi * C.phi_n(0,f) * I3num(n/2, meff))
    denomlist.append(denom)
    numlist.append(num)
    res2 /= denom
    
    mult = 4 * pf * sqrt(pf**2 + meff**2) / (2 * pi**2)
    return (res + res2)*mult

f_f1 = 0.
def f1(n):
    global f_f1
    pf = eos.p_f(n/2)
    mn = C.M[0]
    
    f_f1, = eos.f_eq(np.array([n/2, n/2]), np.array([f_f0]), 1, C)
    f = f_f1
    meff = mn*C.phi_n(0,f)
    
    res = - C.Co * pf**2/(mn**2 * C.eta_o(f) * (pf**2 + meff**2))
    res /= (1. + 2 * C.Co/(mn**2 * C.eta_o(f)) * ( (2./3.) * I1num(n/2, meff) + 
                                                   meff**2 * I4num(n/2, meff)))
    
    mult = 4 * pf * sqrt(pf**2 + meff**2) / (2 * pi**2)
    return res*mult


nlist = np.linspace(0., 4., 100)

klist = []
f0list = []
k1list = []
f1list = []
with open('MatsuiF0.dat', 'r') as f:
    for line in f:
        k, f = line.split()
        klist.append(float(k))
        f0list.append(float(f))

# with open('MatsuiF1.dat', 'r') as f:
#     for line in f:
#         k, f = line.split()
#         k1list.append(float(k))
#         f1list.append(float(f))


n0 = 0.3*wr.n0#2*(1.42 * 197.33 / 135.)**3 /(3 * pi **2)

f_eq0_der = 0.
def f0_der(n):
    global f_eq0_der
    mn = C.M[0]
    pf = eos.p_f(n/2)
    f_eq0_der, = eos.f_eq(np.array([n/2, n/2]), np.array([f_eq0_der]), 1, C)
    meff = mn*C.phi_n(0, f_eq0_der)
    
    print f_eq0_der
    mu = lambda z, f: derivative(lambda x: eos._E(np.array([f, x/2, x/2]), C), z, dx=1e-3)
    
    dmu_df = derivative(lambda x: mu(n, x), f_eq0_der, dx=1e-3)
    dmu_dn = derivative(lambda x: mu(x, f_eq0_der), n, dx=1e-3)
    d2E_df = derivative(lambda z: eos._E(np.array([z, n/2, n/2]), C), f_eq0_der, dx=1e-3, n=2)
    print 'dmu_dn = ', dmu_dn, ' ?= ', C.Co/ mn**2 /C.eta_o(f_eq0_der) + (3 * pi**2)**(2./3) * (n/2)**(-1./3) / sqrt(pf**2 + meff**2) / 6
    print derivative(lambda x: eos._E(np.array([f_eq0_der, x/2, x/2]), C), n, dx=1e-3, n=2)
    res = dmu_dn - (3 * pi**2)**(2./3) * (n/2)**(-1./3) / sqrt(pf**2 + meff**2) / 6 - dmu_df**2 / d2E_df
    
    mult = 4 * pf * sqrt(pf**2 + meff**2) / (2 * pi**2)
    return res*mult

f_eq0_der_NM = 0.
def f0_der_NM(n):
    global f_eq0_der_NM
    mn = C.M[0]
    pf = eos.p_f(n)
    f, = eos.f_eq(np.array([n, 0.]), np.array([f_eq0_der_NM]), 1, C)
    f_eq0_der_NM = f
    
    meff = mn*C.phi_n(0, f)
    
    print f
    mu = lambda z, f: derivative(lambda x: eos._E(np.array([f, x, 0.]), C), z, dx=1e-3)
    
    dmu_df = derivative(lambda x: mu(n, x), f, dx=1e-3)
    dmu_dn = derivative(lambda x: mu(x, f), n, dx=1e-3)
    d2E_df = derivative(lambda z: eos._E(np.array([z, n, 0.]), C), f, dx=1e-3, n=2)
#     print 'dmu_dn = ', dmu_dn, ' ?= ', C.Co/ mn**2 /C.eta_o(f) + (3 * pi**2)**(2./3) * (n/2)**(-1./3) / sqrt(pf**2 + meff**2) / 6
#     print derivative(lambda x: eos._E(np.array([f, x, 0.]), C), n, dx=1e-3, n=2)
    res = dmu_dn - (3 * pi**2)**(2./3) * (n)**(-1./3) / sqrt(pf**2 + meff**2) / 6 - dmu_df**2 / d2E_df
    
    mult = 4 * pf * sqrt(pf**2 + meff**2) / (2 * pi**2)
    return res*mult

n_nm = np.linspace(0., 4., 100)
plt.plot(n_nm, map(f0_der_NM, n_nm))
plt.show()

print f0_der(wr.n0)
print f0(wr.n0)
# exit()

print (135./197.33)*eos.p_f(n0/2)
print 'F_0(n0) = ', f0(n0)
pf0 = eos.p_f(n0/2)
f_0 = 0.
for n in linspace(0, n0, 1000):
    f_0, = eos.f_eq(np.array([n/2, n/2]), np.array([f_0]), 1, C)
#     print f_0
print 'f_0 = ', f_0
print pf0
print C.M[0]
print 'K1 = ', eos.K(n0, C)

print 3 * pf0**2 /sqrt(pf0**2 + (C.M[0]*C.phi_n(0, f_0))**2) * (1 + f0(n0)) * 135.

# wr.solve()
print C.f0
print eos.K(wr.n0, C)
f_0 = C.f0
pf0 = eos.p_f(wr.n0/2)
print 3 * pf0**2 /sqrt(pf0**2 + (C.M[0]*C.phi_n(0, C.f0))**2) * (1 + f0(wr.n0)) * 135.
exit()
denomlist=[]
numlist=[]
fig, ax = plt.subplots(2,1)
ax[0].plot(nlist/wr.n0, map(f0, nlist))
ax[0].plot(nlist/wr.n0, map(f1, nlist))
ax[0].plot(nlist/wr.n0, map(f0_der, nlist))
ax[1].plot(nlist/wr.n0, denomlist)
ax[1].plot(nlist/wr.n0, numlist)
ax[1].plot(nlist/wr.n0, np.array(numlist)/np.array(denomlist))
ax[0].plot(nlist/wr.n0, map(lambda z: f0_der(z) - f0(z), nlist))
# plt.plot(nlist/wr.n0, map(f0, nlist))
# plt.plot(nlist/wr.n0, map(f1, nlist))
ax[0].plot(nlist/wr.n0, [-1. for i in nlist], c='red', ls='--')
ax[0].plot(klist, f0list)
# ax[0].plot(k1list, f1list)

# plt.xlim([0., 5.])
plt.show()

tab0 = np.array([nlist/wr.n0, map(f0, nlist)]).transpose()
table0 = tabulate(tab0, tablefmt='plain')

tab1 = np.array([nlist/wr.n0, map(f1, nlist)]).transpose()
table1 = tabulate(tab1, tablefmt='plain')

# with open(fname_0, 'w') as f:
#     f.write(table0)
# 
# with open(fname_1, 'w') as f:
#     f.write(table1)    



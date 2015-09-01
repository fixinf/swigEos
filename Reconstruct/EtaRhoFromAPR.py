import eosWrap as eos
import numpy as np
import matplotlib.pyplot as plt
from Wrapper import Wrapper, FromEos
from scipy import integrate
from pylab import pi, sqrt
from scipy.interpolate import interp1d, interpolate
from scipy.misc.common import derivative
from tabulate import tabulate
from scipy.interpolate.fitpack2 import UnivariateSpline
C = eos.KVOR()
IC = eos.InterpolatedScalings()

wr = Wrapper(C)
iwr = Wrapper(IC)

IC.Cs = C.Cs
IC.Co = C.Co
IC.Cr = C.Cr
IC.b = C.b
IC.c = C.c

EC = FromEos()
EC.setConstants(Co=1*C.Co, Cr=IC.Cr)
EC.setParams()
# EC.setScaling(eta_o=lambda z: ((1 + 0.65 * z) / (1 + 0.65 * 0.15), 
#               eta_r=lambda z: 1.)

cAPR = 1.0/0.16
nAPR=cAPR*np.array([0.0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.32, 0.4, 0.48, 0.56, 0.64, 0.8, 0.96])
APR = np.array([0.0, -6.48, -12.13, -15.04, -16.00, -15.09, -12.88, -5.03, 2.13, 15.46, 
             34.39, 58.35, 121.25, 204.02])

APR_N = np.array([0., 4.45, 6.45 , 9.65, 13.29, 17.94, 22.92, 27.49, 38.82, 54.95, 75.13, 99.75, 127.58, 205.34, 305.87])
nAPR_N = cAPR* np.array([0, 0.02, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24,  0.32,  0.4, 0.48, 0.56, 0.64, 0.8, 0.96]) 

nAPR = nAPR_N
APR = APR_N

iAprE = interp1d(nAPR, APR, kind='cubic')
# iAprE = UnivariateSpline(nAPR, APR)

aprEps = nAPR*wr.n0*APR/wr.m_pi + C.M[0]*nAPR*wr.n0


i_apr_eps = interp1d(nAPR*wr.n0, aprEps, kind='cubic')
# i_apr_eps = UnivariateSpline(nAPR*wr.n0, aprEps)

i_apr_eps2 = lambda z: z*iAprE(z/wr.n0)/wr.m_pi + C.M[0]*z
i_apr_P = lambda z: z*derivative(i_apr_eps, z, dx=1e-5, order=7) - i_apr_eps(z)
i_apr_P2 = lambda z: z*derivative(i_apr_eps2, z, dx=1e-5, order=7) - i_apr_eps2(z)

n = np.linspace(0, 1.*wr.n0, 100)

print(i_apr_P(wr.n0), i_apr_P2(wr.n0))

# plt.plot(n, i_apr_eps(n), n, i_apr_eps2(n))
plt.plot(n[1:], i_apr_P(n[1:]), n[1:], i_apr_P2(n[1:]))
plt.show()

i_apr_eps = i_apr_eps2
i_apr_P = i_apr_P2

plt.plot(n, [135*(i_apr_eps(z)/z - C.M[0]) for z in n], n, 135*(wr.Eneutr(n)/n - C.M[0]))
# plt.plot(n, iAprE(n/wr.n0), n, 135*(wr.Eneutr(n)/n - C.M[0]))
plt.show()

E, fn0 = wr.Eneutr(n, ret_f=1)
P = wr.P_N(n) / wr.const / wr.m_pi**4

E = i_apr_eps(n)
P = i_apr_P(n[1:])
P = np.insert(P, 0, 0.)


# EC.setFs(fs0, n)
# plt.plot(n, E, n, map(EC.ES, n))
# plt.show()

Fn = EC.FnFromEos(E, P, n)

plt.plot(n, fn0, n, Fn)
plt.show()

Ulist = []

for i, _f in enumerate(Fn):
    res = - IC.M[0]**4 * _f**2 / 2 / IC.Cs
    res -= eos.kineticInt(n[i], C.M[0]*(1 -_f), _f)
    res += integrate.quad(lambda z: sqrt(EC.pf(z)**2 +
             EC.mn**2*(1-EC.Fn(z))**2), 0., n[i], epsrel=1e-18)[0]
    Ulist.append(res)
    
Ulist = np.array(Ulist)
UKV = np.array(list(map(C.U, fn0)))
lines=plt.plot(Fn, Ulist, fn0, list(map(C.U, fn0)))
plt.ylabel('U', fontsize=24)
plt.xlabel('f', fontsize=24)
plt.legend(lines, ['Reconstructed', 'KVOR'], fontsize=24, loc=0)
plt.show()

tab = np.array([Fn, Ulist]).transpose()
table = tabulate(tab, ['f', 'U [m_pi^4]'], tablefmt='plain')
with open('U_APR_Co=%.2f.dat'%EC.Co, 'w') as f:
    f.write(table)


# Ulist.insert(0, 0.)

nf = interp1d(Fn, n)

eta_o = IC.Cr * n ** 2 /(8 * C.M[0]**2) / (Ulist - UKV + EC.Cr * n ** 2 /(2 * C.M[0]**2))
plt.plot(Fn, eta_o)
plt.show()

# IC.set_eta_s(np.insert(Fs, 0, -1e-3), np.insert(Ulist, 0, 0.))
# IC.set_eta_o(np.linspace(-1e-3, 1, 100), np.linspace(1., 1+1e-9, 100))

# eta_o = np.ones(Fs.shape)

# eta_o = map(C.eta_o, Fs)
print(max(Fs))
fig, ax = plt.subplots(2,1)
ax[0].plot(n, eta_o, n, list(map(C.eta_o, fs0)))
ax[1].plot(n, [135*(i_apr_eps(z)/z - C.M[0]) for z in n], n, 135*(wr.Esymm(n)/n - C.M[0]))
plt.show()

# eta_o[10] += 0.01
print(eta_o)
# exit()
eta_o[0] = 0.55
IC.set_eta_s(np.insert(Fs, 0, -1e-3), np.insert(UKV, 0, 0.))
IC.set_eta_o(np.insert(Fs, 0, -1e-3), np.insert(eta_o, 0, 0.5))
# IC.set_eta_o(Fs[1:], eta_o[1:])

IC.fmax = np.max(Fs)
IC.Co = EC.Co

n_i = n[5:-1]

iE, iF = iwr.Esymm(n_i, ret_f=1)
print(iF)
plt.plot(n_i, iF, n, Fs)
plt.show()
plt.plot(Fs, Ulist, iF, list(map(IC.U, iF)), iF, list(map(C.U, iF)))
plt.show()
plt.plot(iF, list(map(IC.eta_o, iF)), iF, list(map(C.eta_o, iF)))
plt.show()

eta_tab = np.array([iF, list(map(IC.eta_o, iF)), list(map(C.eta_o, iF))]).transpose()
f_tab = np.array([n_i, iF, n, Fs]).transpose()

plt.plot(n_i, iE, n, E)
plt.show()
IC.f0 = 0.195
# print iwr.K()
Ebind = (iE/n_i - C.M[0])*wr.m_pi
EbindOrig=(E / n - C.M[0])*wr.m_pi
EKV = (wr.Esymm(n) / n - C.M[0])*wr.m_pi



e_tab = np.array([n_i/wr.n0, Ebind, n/wr.n0, EbindOrig, EKV]).transpose()
plt.plot(n_i/wr.n0, Ebind, n/wr.n0, EbindOrig, n/wr.n0, EKV)
plt.show()

eta_table = tabulate(eta_tab,['f','Interpolated','KVOR'], tablefmt='plain')
f_table = tabulate(f_tab, ['n_inter','iF', 'n', 'fsAPR'], tablefmt='plain')
e_table = tabulate(e_tab, ['n_i', 'ebind', 'n', 'ebind_orig','ebind KVOR'], tablefmt='plain')

with open('eta.dat', 'w') as f:
    f.write(eta_table)
    
with open('f.dat', 'w') as f:
    f.write(f_table)
    
with open('e.dat', 'w') as f:
    f.write(e_table)
    

import eosWrap as eos
import numpy as np
import matplotlib.pyplot as plt
from Wrapper import Wrapper, FromEos
from scipy import integrate
from pylab import pi, sqrt
from scipy.interpolate import interp1d
from scipy.misc.common import derivative
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
EC.setConstants(Co=1.3*IC.Co, Cr=IC.Cr)
EC.setParams()

cAPR = 1.0/0.16
nAPR=cAPR*np.array([-1e-2, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.32, 0.4, 0.48, 0.56, 0.64, 0.8, 0.96])
APR = np.array([0.0, -6.48, -12.13, -15.04, -16.00, -15.09, -12.88, -5.03, 2.13, 15.46, 
             34.39, 58.35, 121.25, 204.02])

aprEps = nAPR*wr.n0*APR/wr.m_pi + C.M[0]*nAPR*wr.n0


i_apr_eps = interp1d(nAPR*wr.n0, aprEps, kind='cubic')
i_apr_P = lambda z: z*derivative(i_apr_eps, z, dx=1e-3) - i_apr_eps(z)

n = np.linspace(0., 5.5*wr.n0, 100) 
E, fs0 = wr.Esymm(n, ret_f=1)
P = wr.Psymm(n) / wr.const / wr.m_pi**4

# E = i_apr_eps(n)
# P = i_apr_P(n)


# EC.setFs(fs0, n)
# plt.plot(n, E, n, map(EC.ES, n))
# plt.show()


Fs = EC.FsFromEos(E, P, n)
plt.plot(n, fs0, n, Fs)
plt.show()

Ulist = []

for i, _f in enumerate(Fs):
    res = - IC.M[0]**4 * _f**2 / 2 / IC.Cs
    res -= 2*eos.kineticInt(n[i]/2, C.M[0]*(1 -_f), _f)
    res += integrate.quad(lambda z: sqrt(EC.pf(z/2)**2 +
             EC.mn**2*(1-EC.Fs(z))**2), 0., n[i], epsrel=1e-18)[0]
    Ulist.append(res)
Ulist = np.array(Ulist)

# Ulist.insert(0, 0.)

lines=plt.plot(Fs, Ulist, fs0, map(C.U, fs0))
plt.legend(lines, ['Reconstructed', 'KVOR'])
plt.show()

IC.set_eta_s(np.insert(Fs, 0, -1e-3), np.insert(Ulist, 0, 0.))
IC.fmax = np.max(Fs)
IC.Co = EC.Co
iE, iF = iwr.Esymm(n[1:-1], ret_f=1)
plt.plot(n[1:-1], iF, n, Fs)
plt.show()
plt.plot(Fs, Ulist, iF, map(IC.U, iF))
plt.show()
plt.plot(n[1:-1], iE, n, E)
plt.show()

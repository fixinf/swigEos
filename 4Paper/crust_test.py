from os.path import join

from scipy.interpolate._monotone import Akima1DInterpolator
from scipy.interpolate.fitpack2 import InterpolatedUnivariateSpline, UnivariateSpline
from scipy.interpolate.interpolate import interp1d

import Models2

__author__ = 'const'
import matplotlib
matplotlib.use("QT5Agg")
from matplotlib import pyplot as plt
from RhoCond.rho_wrap import *
import numpy as np

wr = Models2.KVOR()
# wr.nucl.stars_crust(show=1, inter='cubic')
# exit()
# wr.dumpAll(hyper=0)
data_KV = np.loadtxt(join(wr.foldername, wr.filenames['eos']), skiprows=1)

NKV = data_KV[:, 0] * wr.n0
EKV = data_KV[:, 4]
PKV = data_KV[:, 5]

name = 'KVOR'

NN, EN, PN = np.loadtxt(name + '_4EE.dat').transpose()

data2 = np.loadtxt('pe_'+name+'.dat')

Ec, Pc, Nc = np.loadtxt('crust_export.dat').transpose()

# plt.plot(Nc, Ec/Nc/wr.n0)
# plt.plot(NN/wr.n0, EN/NN)
# plt.plot(NKV/wr.n0, EKV/NKV)
# plt.show()


n_cr = 0.3
n_eos = 0.9

n_low = np.array([n for n in Nc if n < n_cr])
n_high = np.array([n for n in NKV if n > n_eos*wr.n0])

print(n_low.shape, n_high.shape)

n_inter = np.concatenate((n_low*wr.n0, n_high))
e1 = Ec[:n_low.shape[0]]/n_low/wr.n0
e2 = EKV[EKV.shape[0] - n_high.shape[0]:]/n_high

e_inter = np.concatenate((e1, e2))

print(n_inter)

nrange = np.linspace(0., 3*wr.n0, 500)

# iE = interp1d(n_inter, np.nan_to_num(e_inter), kind='cubic', bounds_error=0, fill_value=0.)
iE = UnivariateSpline(n_inter, np.nan_to_num(e_inter), s=1e-3, k=3)

ielist = iE(nrange)
print(ielist)

plt.plot(nrange/wr.n0, ielist)
plt.plot(Nc, Ec/Nc/wr.n0)
# plt.plot(NN/wr.n0, EN/NN)
plt.plot(NKV/wr.n0, EKV/NKV)
plt.xlim([0., 2.])
plt.ylim([6.8, 7.4])

plt.show()


deriv_P = nrange * nrange * np.array([derivative(iE, z, dx=1e-6) for z in nrange])
plt.plot(nrange/wr.n0, deriv_P)
plt.plot(Nc, Pc)
plt.plot(NKV/wr.n0, PKV)
plt.ylim([0., 0.2])
plt.xlim([0., 2])
plt.show()

i_eps = ielist*nrange
plt.plot(i_eps, deriv_P)
plt.plot(EN, PN)
plt.show()


# np.savetxt(name + "_newcrust.dat", np.array([nrange/wr.n0, ielist, i_eps, deriv_P]).transpose(),
#            fmt='%.6f')

# plt.plot(data2[:, 0], data2[:, 1], 'bs')
# plt.xlim([0, 0.8])
plt.show()


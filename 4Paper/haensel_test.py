from os.path import join

from scipy.integrate.quadpack import quad
from scipy.interpolate.interpolate import interp1d
from scipy.optimize._root import root

import Models2

__author__ = 'const'
import matplotlib
matplotlib.use("QT5Agg")
from matplotlib import  pyplot as plt
from RhoCond.rho_wrap import *
import numpy as np


wr = Models2.myMod()
# wr.nucl.stars_crust(show=1, inter='cubic')
# exit()
# wr.dumpAll(hyper=0)
data_KV = np.loadtxt(join(wr.foldername, wr.filenames['eos']), skiprows=1)

NKV = data_KV[:, 0] * wr.n0
EKV = data_KV[:, 4]
PKV = data_KV[:, 5]

name = 'MKVOR'
print(wr.n0)
# exit()
NN, EN, PN = np.loadtxt(name + '_4EE.dat').transpose()

data2 = np.loadtxt('pe_'+name+'.dat')

Ec, Pc, Nc = np.loadtxt('crust_export.dat').transpose()
Nc *= wr.n0


#####Working with EN, PN (no proper python interpolation for P(E) for MKVOR)
kind = 'linear'

iPEcrust = interp1d(Ec, Pc, kind=kind)
iENcrust = interp1d(Nc, Ec, kind=kind)
iNEcrust = interp1d(Ec, Nc, kind=kind)
iPNcrust = interp1d(Nc, Pc, kind=kind)

iPEkv = interp1d(EKV, PKV, kind=kind)
iNEkv = interp1d(EKV, NKV, kind=kind)

################ stupid joining


iPE = interp1d(EN, PN, kind=kind)
iEN = interp1d(NN, EN, kind=kind)
iNE = interp1d(EN, NN, kind=kind)
iPN = interp1d(NN, PN, kind=kind)


n_s = 0.7 * wr.n0
E_s = iEN(n_s)

E_join = root(lambda z: iPEkv(z) - iPEcrust(z), iENcrust(wr.n0))['x']

erange = np.linspace(E_s, 0, 200, endpoint=1)
P_stupid = np.array([iPEcrust(z) for z in erange[::-1] if z < E_join]  +
                    [iPEkv(z) for z in erange[::-1]  if z > E_join ])[::-1]
print(P_stupid)
print(E_s)
# plt.plot(erange, iPE(erange))
# plt.plot(erange, P_stupid)
# plt.plot(erange, iPEkv(erange))
# plt.plot(erange, iPEcrust(erange))
# plt.legend()
# plt.show()

iPE = interp1d(erange, P_stupid)

#####Stupid raw sum

func = lambda z: 1. / (iPE(z) + z)

# def func(z):
#     print(z)
#     return 1. / (iPE(z) + z)

print(n_s, E_s)

nrange = np.array([quad(func, E_s, z)[0] for z in erange])
nrange_test = np.array([iNE(z) for z in erange])

nrange = n_s * np.exp(nrange)
iNEnew = interp1d(erange, nrange)
# nrange = np.insert(nrange, nrange.shape[0], 0.)
print(nrange)
# exit()
Nnew = iNEnew(erange)
Neos = iNEkv(erange)
Ncrust = iNEcrust(erange)

Pnew = iPE(erange)
Peos = iPEkv(erange)
Pcrust = iPEcrust(erange)

iPnew = interp1d(Nnew, Pnew, bounds_error=0, fill_value=0.)
iEnew = interp1d(Nnew, erange, bounds_error=0, fill_value=0.)
# plt.plot(Nnew, Pnew)
# plt.plot(NN, PN)
# plt.show()

# plt.plot(Nnew, erange)
# plt.plot(NN, EN)
# plt.show()

nmin = 0.#min(nrange)

pdiff = quad(lambda z: abs(iPN(z) - iPnew(z)), nmin, n_s)[0] / (n_s /
                    wr.n0 * 0.16) * wr.mpi4_2_mevfm3


ediff = quad(lambda z: abs(iEN(z) - iEnew(z)), nmin, n_s)[0] / (n_s /
                    wr.n0 * 0.16) * wr.mpi4_2_mevfm3

epsdiff = quad(lambda z: abs(iEN(z)/z - iEnew(z)/z), nmin, n_s)[0] / (n_s /
                    wr.n0) * wr.m_pi


pdiff_rel = quad(lambda z: abs(iPN(z) - iPnew(z))/iPN(z), nmin, n_s)[0] / (n_s /
                    wr.n0)

ediff_rel = quad(lambda z: abs(iEN(z) - iEnew(z))/iEnew(z), nmin, n_s)[0] / (n_s /
                    wr.n0)


epsdiff_rel = quad(lambda z: abs(iEN(z)/z - iEnew(z)/z) / (iEnew(z)/z), nmin, n_s)[0] / (n_s /
                    wr.n0)

print('p', pdiff)
print('e', ediff)
print('eps', epsdiff)
print('p_rel', pdiff_rel)
print('e_rel', ediff_rel)
print('eps_rel', epsdiff_rel)



exit()

out = np.array([erange * wr.mpi4_2_mevfm3,
                Pnew * wr.mpi4_2_mevfm3,
                Nnew / wr.n0,
                Peos * wr.mpi4_2_mevfm3,
                Neos / wr.n0,
                Pcrust * wr.mpi4_2_mevfm3,
                Ncrust / wr.n0,
                erange/Nnew * wr.m_pi,
                erange/Neos * wr.m_pi,
                erange/Ncrust * wr.m_pi]).transpose()


np.savetxt('outKVOR_rude.dat', out, fmt='%.8f', header='# E | P_new | N_new | P_kv | N_kv | P_bps | N_bps |' +
           ' eps_new | eps_kv | eps.crust')



plt.plot(nrange, label='result')
plt.plot(nrange_test, label='test')
plt.legend()
plt.show()

plt.plot(nrange_test[:] / wr.n0, (nrange - nrange_test)/nrange_test, label='\eps')
plt.legend()
plt.show()

np.savetxt(name+'_n_out.dat', np.array([nrange_test/wr.n0, nrange/wr.n0, erange, iPE(erange)]).transpose(),
           fmt = '%.6f')

plt.plot(nrange_test, erange/nrange)
plt.plot(nrange_test, erange/nrange_test)
plt.plot(NKV, EKV/NKV)
plt.show()
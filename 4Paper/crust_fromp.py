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

nrange = np.linspace(0, 2*wr.n0, 100)

epsnew = wr.C.M[0] + np.array([quad(lambda z: iPN(z)/z**2, 0., x)[0] for x in nrange])
Enew = nrange * epsnew
print(epsnew.shape, NN.shape)
plt.plot(nrange, Enew)
plt.plot(NN, EN)
plt.show()


out = np.array([
    nrange/wr.n0,
    iEN(nrange)*wr.mpi4_2_mevfm3,
    iPN(nrange)*wr.mpi4_2_mevfm3,
    Enew*wr.mpi4_2_mevfm3
]).transpose()

np.savetxt(name+'_from_p.dat', out, fmt = '%.8f', header='n [n_0] | E_old [mev/fm3] '
                                                         '| P_old [..] | E_new [..]')
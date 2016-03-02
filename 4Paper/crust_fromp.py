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


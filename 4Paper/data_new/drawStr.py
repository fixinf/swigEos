__author__ = 'const'
import Models2
from scipy import interpolate
from tabulate import tabulate
from os.path import join
import numpy as np
from numpy import array as arr

m = Models2.KVOR()

def getProfile(n_c, m=Models2.KVOR().hyper):
    m.dumpMassesCrust(nmin=n_c, nmax=n_c, npoints=1, write=False)
    lastN = m.dr.getLastN(m.dr.nSize)[:-1]
    lastR = m.dr.getLastR(m.dr.nSize)[:-1]
    lastM = m.dr.getLastM(m.dr.nSize)[:-1]
    grav_mult = []
    for i, r in enumerate(lastR):
        grav_mult.append(1. / np.sqrt(1 - 2 * 1.4677 * lastM[i] / r))
    ###This grav_mult differs by a lack of r^2 from that in dumpMassesCrust
    grav_mult = np.array(grav_mult)

    dx = lastR[1] - lastR[0]
    conc = m.concentrations()
    inter_hyp = [interpolate.interp1d(m.nrange, conc[:, i])
                 for i in range(1, m.n_baryon)]
    return lastN, lastR, grav_mult, arr(map(lambda z: [f(z) for f in inter_hyp], lastN)).transpose().tolist()

def writeProfile(m, n_c, fname):
    n, r, gm, part = getProfile(n_c, m)
    with open(join(m.foldername, fname), 'w') as f:
        f.write(tabulate(arr([n/m.n0, r] + part + [gm]).transpose(), ['n/n0', 'r[km]']
                         + m.part_names[1:] + ['gamma'], tablefmt='plain'))

writeProfile(m.hyper, 6.18 * m.n0, 'mmax_hyper.dat')
writeProfile(m.hyper_phi, 7.24 * m.n0, 'mmax_hyper_phi.dat')
writeProfile(m.hyper_phi_sigma, 6.71532*m.n0, 'mmax_hyper_phi_sigma.dat')






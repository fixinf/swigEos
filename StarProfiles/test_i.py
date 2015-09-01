from math import pi
from scipy.integrate.odepack import odeint
from scipy.interpolate.interpolate import interp1d

__author__ = 'const'
import numpy as np
import Models2
from matplotlib import pyplot as plt

def func_i(y, t, E, P, M):
    print(t)
    mpi4 = 5.7207e-5
    res = mpi4 * (8 * pi/3) * t**4 * (1 - 2 * y / t**3) * (1 - y / 2 / t**3) * (E(t) +
            P(t)) / (1 - 2 * 1.4766 * M(t) / t)
    # print(res)
    return res

def get_i(m):
    res =  m.dumpMasses(nmin=1.5, nmax=1.5, npoints=1, write=0)
    nsize = m.dr.nSize
    R = m.dr.getLastR(nsize)[:-1]
    P = m.dr.getLastP(nsize)[:-1]
    E = m.dr.getLastE(nsize)[:-1]
    M = m.dr.getLastM(nsize)[:-1]
    # plt.plot(R, M)
    # plt.show()
    fill = 0.
    kind = 'cubic'
    args = (interp1d(R, E, bounds_error=0, fill_value=E[-1], kind=kind),
            interp1d(R, P, bounds_error=0, fill_value=P[-1], kind=kind),
            interp1d(R, M, bounds_error=0, fill_value=M[-1], kind=kind)
    )

    rRange = np.linspace(R[1], R[-1], 100)
    print(res)
    res = odeint(func_i, [0], rRange, args=args, )
    print(res[-1] / 1.4766 * 2e33 * 1e10)

m = Models2.KVOR().nucl
get_i(m)
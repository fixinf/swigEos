from os.path import join
from scipy.interpolate.interpolate import interp1d
from scipy.optimize.zeros import bisect

__author__ = 'const'
import Models2
from matplotlib import pyplot as plt
import numpy as np
plist = [[0.8, 7.5] ,
        [1.2, 7.5],
         [1.5, 7.5],
         [1.4, 9.]
         ]

plist = [[0.,0.], [0.8, 3.6], [0.8, 6]]

def mass_f(n, m_target, m=Models2.KVOR()):
    nrange, mg ,r, mb1, mb2, = m.stars_crust(nmin=n, nmax=n, npoints=1)
    print(n/m.n0, mg[0] - m_target)
    return mg[0] - m_target

for p in plist:
    m = Models2.myMod_L(*p)
    n_c = bisect(lambda z : mass_f(z, 1.25, m), 2*m.n0, 3*m.n0, xtol=1e-4)
    print(n_c, mass_f(n_c, 1.25, m))
    m.nucl.dumpProfiles(n_c)

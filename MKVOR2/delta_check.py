from scipy.interpolate.interpolate import interp1d
from scipy.optimize.minpack import leastsq

import Models2
import Models
from matplotlib import pyplot as plt
import numpy as np

# wr = Models2.myModExpOmega(0.68)
wr = Models2.MKValpha03(0.68)
wr.dumpAll(hyper=0)
# wr = Models2.MKVOR_d()
n =  2.662*wr.n0
nd = 0.015083*wr.n0
f = 0.642106
m = wr.delta_sym
S,V = wr.dumpPotentials()
i = 8

xo=1.
U = -125.

iS = interp1d(m.nrange/m.n0, S)
iV = interp1d(m.nrange/m.n0, V)

xs_d = (U - m.C.X_o[i]*iV(1.))/iS(1.)
# exit()

wr.setDeltaConst(np.array([xs_d for i in range(4)]),
             np.array([xo for i in range(4)]),
             np.array([1., 1., 1., 1.]),
             's = %.2f o = %.2f U = %.0f' % (xs_d, xo, U))

wr.delta_sym.checkEq(n, nd, f)
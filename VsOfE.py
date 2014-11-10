import eosWrap as eos
import matplotlib
from scipy.misc.common import derivative
from tabulate import tabulate
matplotlib.use('QT4Agg')
from Wrapper import Wrapper
from scipy.integrate.quadpack import quad
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import Models

C = Models.KVOR()

wr = Wrapper(C)

wr.reset(nmin=0., npoints=500)

E, P, N = wr.EPN()
E *= wr.const

dE = np.diff(E)
dP = np.diff(P)*wr.const

E2 = E[np.argmin(abs(N - [2*wr.n0 for i in N]))]
plt.figure(figsize=(6,6))
line = plt.plot(E[1:], dP/dE)
plt.xlim([0., 800.])
plt.gca().axvline(E2)
plt.xlabel(r'$E \, [MeV/fm^3]$', fontsize=18)
plt.ylabel(r'$v_S^2$',fontsize=18)
plt.legend(line, ['KVOR'], loc=0,fontsize=18)
plt.show()
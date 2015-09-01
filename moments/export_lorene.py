from scipy.integrate.quadpack import quad
from scipy.interpolate.interpolate import interp1d

__author__ = 'const'
import Models2
import numpy as np

m = Models2.KVOR()
npoints = 100


mevfm3 = 1.7827e12 # g/cm^3
mevfm3_dyn = 1.6022e33

nrange = np.linspace(0, 4, npoints)
E, P, n = m.neutr.EPN(nrange=nrange)
# iE = interp1d(P, E)

E *= m.mpi4_2_mevfm3
E *= mevfm3
P *= mevfm3_dyn
n = n / m.n0 * 0.16e39


# H = np.array(list((map(lambda p: quad(lambda z: c**2/(iE(z) + z), 1e-3, p)[0], P))))

f = open(m.Ctype.__name__ + '_lorene.dat', 'w')
# f.write(str(npoints) + '\n')
for i in range(npoints):
    f.write('%i %e %e %e \n' % (i+1, n[i], E[i], P[i]))

f.close()

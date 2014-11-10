import matplotlib
matplotlib.use('QT4Agg')
from pylab import *
import numpy as np
import matplotlib.pyplot as plt


a = []
with open('BSk19_a.dat', 'r') as f:
    for line in f:
        _a, = line.split()
        a.append(float(_a))
        
def zeta_exp(a1, a2, xi):
    return 1./(exp(a1*(a2-xi)) + 1)

def zeta(xi, a):
    res = (a[0] + a[1]*xi + a[2]*xi**3)*zeta_exp(a[4], -a[5], -xi)/(1 + a[3]*xi)
    res += (a[6] + a[7]*xi) * zeta_exp(a[8], a[5], xi)
    res += (a[9] + a[10]*xi) * zeta_exp(a[11], a[12], xi)
    res += (a[13] + a[14]*xi) * zeta_exp(a[15], a[16], xi)
    res += a[17]/(1 + (a[18]*(xi - a[19]))**2)
    res += a[20]/(1 + (a[21]*(xi - a[22]))**2)
    return res
print size(a)
xilist = np.linspace(0., 10., 100)
plt.plot(exp(xilist), exp(map(lambda z: zeta(z, a), xilist)))
plt.show()
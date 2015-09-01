__author__ = 'const'

from matplotlib import pyplot as plt
import Models2
import numpy as np
import numpy.polynomial.chebyshev as cheb
import numpy.polynomial.hermite as herm

x = np.linspace(0, 6, 100)

m = Models2.KVOR()

E, P, n = m.nucl.EPN()
P /= m.nucl.mpi4_2_mevfm3
vs = m.nucl.dumpVs(write=False)
plt.plot(n/m.n0, vs)
plt.show()
gamma = np.nan_to_num((1 + E/P) * vs)

plt.plot(n/m.n0, gamma)
plt.show()

# n = np.linspace(0, 3, 100)
# gamma = np.sin(n)

mean = np.trapz(gamma, n) / (n[-1] - n[0])
print (mean)
gamma = gamma - mean
# exit()
ch = cheb.chebfit(n, gamma, 5)
print(ch)
plt.plot(n/m.n0, gamma)
plt.plot(n/m.n0, cheb.chebval(n, ch))
# plt.show()

ht = herm.hermfit(n, gamma, 5)
print(ht)
# plt.plot(n/m.n0, gamma)
plt.plot(n/m.n0, herm.hermval(n, ch))
plt.show()

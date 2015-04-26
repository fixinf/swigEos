__author__ = 'const'
import Models2
from matplotlib import pyplot as plt
import numpy as np
m = Models2.KVOR()

E, P, n = m.nucl.EPN()
Eparts = m.nucl.Eparts

plt.plot(n, E)
plt.plot(n, Eparts)
print Eparts
# plt.plot(n, np.sum(Eparts, axis=1).transpose())
plt.show()
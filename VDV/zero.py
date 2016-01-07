from os.path import join

__author__ = 'const'
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import Models2
wr = Models2.KVOR()
# wr.dumpEos()
# data = np.loadtxt(join(wr.foldername, wr.filenames['eos']),skiprows=1)
E, P ,n = wr.sym.EPN()
# n = data[:,0]
# P = data[:,2]
V = 0.16 / n
plt.plot(V, P)
plt.ylim([-10,10])
# plt.xlim([0,150.])
plt.show()
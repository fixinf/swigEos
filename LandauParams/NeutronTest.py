import matplotlib
matplotlib.use("QT4Agg")
import Models
from Wrapper import Wrapper
from matplotlib import pyplot as plt
import numpy as np

C = Models.myMod()
wr = Wrapper(C)

data = np.loadtxt('/home/const/GrabbedFigures/MatsuiNeutron/MatsuiNeutronConverted.dat', delimiter=' ', skiprows=1)
n = data[:, 0]*wr.n0
n = np.linspace(0., 3., 400)
# f0 = wr.f0_nm(n)
# plt.plot(data[:,0], data[:,1], n/wr.n0, f0)

f1 = wr.f1_nm(n)
plt.plot(data[:,0], data[:,2], n/wr.n0, f1)
plt.show()
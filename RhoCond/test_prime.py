import Models2
from matplotlib import pyplot as plt
import numpy as np
wr = Models2.KVOR()


m = wr.rcp_nucl


# frange = np.linspace(0, 1., 100)
# plt.plot(frange, [m.C.chi_prime(f) for f in frange])
# plt.show()
# exit()
m.dumpEos()
data = np.loadtxt('/home/const/Dropbox/GrabbedFigures/KV2004/np_chi.csv',
    delimiter=',', skiprows=1)
plt.plot(data[:, 0], data[:, 1])
plt.plot(m.nrange/m.n0, m.concentrations())
plt.plot(m.nrange/m.n0, m.nc/m.nrange)
plt.show()
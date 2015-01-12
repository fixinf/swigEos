import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join


model = Models.myMod

f0 = 0.195
K = 240.

folderName = join('/home/const/Dropbox/Documents/For DN/Very final/data',
                  model.__name__)

if model.__name__ == 'KVORLowerK':
    folderName += '%.2f_%.0f'%(f0, K)

if model.__name__ == 'myModLowerK':
    folderName += '%.0f'%(K)




if model.__name__ == 'KVORLowerK':
    args = (f0, K)
elif model.__name__ == 'myModLowerK':
    args = (K,)
else:
    args=()
        
C = model(*args)

wr = Wrapper(C)
C.Hyper=0
wr.dumpStarProfile(5.67, folderName)
exit()
# wr.reset()
# wr.setDriver()
n_star = 5.25*wr.n0

n, m, r, mg1, mg2 = wr.stars_crust(nmin=wr.n0, nmax=n_star+0.01, npoints=2)

E = wr.dr.getLastE(wr.dr.nSize)
N = wr.dr.getLastN(wr.dr.nSize)
P = wr.dr.getLastP(wr.dr.nSize)
R = wr.dr.getLastR(wr.dr.nSize)
# plt.plot(N/N[0], E)
# plt.plot(N/N[0], P)
# plt.plot(N/N[0], N/N[0], ls='--')

plt.plot(R, E)
plt.plot(R, P)
ax = plt.gca()
# ax.invert_xaxis()
plt.show()
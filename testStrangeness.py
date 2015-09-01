import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join

model = Models.KVOR
C = model()
folderName = join('/home/const/Dropbox/Documents/For DN/Very final/data',
                      model.__name__)
# C = Models.KVOR_cut_03()

wr = Wrapper(C)
n, m, r, mg1, mg2, frac = wr.dumpMassesCrustHyper(folderName, ret_frac=1, fasthyp=1)
i = np.argmax(m)
print('hyp_frac(max) = ', frac[i],'m_max = ', m[i])
plt.plot(n/wr.n0, frac)
plt.show()
plt.plot(frac, m)
plt.show()

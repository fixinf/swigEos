import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join


models = [Models.Cubero, Models.Cubero_cut03, Models.Cubero_cut04, Models.Cubero_cut02]
# models = [Models.Cubero_cut02]
# models = [Models.Cubero]

DUs = []

for model in models:
    args=()

    C = model(*args)
    C.set_xs(np.array([0., 0., -28., 30., 30., 30., -15., -15.]))
    wr = Wrapper(C)
#     exit()

#     n = np.linspace(0, 32*wr.n0, 200)
#     E, f = wr.Esymm(n, ret_f=1)
#     plt.plot(n/wr.n0, f)
#     plt.show()
#     wr.testHyperBind()

    folderName = join('/home/const/Dropbox/Documents/For DN/Very final/data',
                      model.__name__)

    print folderName
    C.Hyper=0
#     wr.dumpMasses(folderName)

    DUs.append([model.__name__, wr.eval_DU(crust=1)])
    wr.dumpAll(folderName, nmax=9.)
#     wr.dumpMassesCrust(folderName, show=0, neutron=0, ncut_crust = .6, ncut_eos = .9, inter='cubic')

print DUs

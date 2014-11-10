import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join


eta_s_list = []
models = [Models.KVOR_cut_smooth]
for model in models:
    C = model()
    C.set_xs(np.array([0., 0., -28., 30., 30., 30., -15., -15.]))
    wr = Wrapper(C)
    
#     wr.testHyperBind()
    
    folderName = join('/home/const/Dropbox/Documents/For DN/Very final/data',
                      model.__name__)
#     wr.dumpEos(folderName)
#     exit()
    print folderName
    wr.dumpAll(folderName)
    wr.reset(npoints=100)
    eta_s_list.append(map(lambda z: C.eta_s(z), wr.rho[:, 0]))
    #wr.dumpEos(folderName)
    
eta_s_list = np.array(eta_s_list).transpose()
plt.plot(wr.n/wr.n0, eta_s_list)
plt.show()


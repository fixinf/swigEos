import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join
from time import sleep


# models = [Models.KVOR, Models.KVOR_cut_mod, Models.myMod]
# models = [Models.myMod3]
models = [Models.myModDiffL]
gamma = 9.5
beta = 1.2
for model in models:
    C = model(beta=beta, gamma=gamma)
    C.set_xs(np.array([0., 0., -28., 30., 30., 30., -15., -15.]))
    wr = Wrapper(C)
    
    print wr.L()
    sleep(1)
#     wr.testHyperBind()
    
    folderName = join('/home/const/Dropbox/Documents/For DN/Very final/data',
                      model.__name__+
                      '%.2f_%.2f'%(beta, gamma))
#     wr.dumpEos(folderName)
#     exit()
    print folderName
    
#     wr.dumpMasses(folderName)
#     wr.dumpEos(folderName)
#     wr.dumpAll(folderName, folderName+'.zip')
#     wr.dumpScalings(folderName)
    wr.dumpMassesCrust(ncut_crust=0.6, ncut_eos=0.8, folderName=folderName)
#     wr.dumpLandauParams(folderName)
    print wr.L()

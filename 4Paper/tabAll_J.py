import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join


models = [Models.KVOR, Models.KVOR_cut_mod, Models.myMod]
# models = [Models.myMod3]
# models = [Models.KVOR_cut_mod]
# models = [Models.KVOR_cut_02]
# models = [Models.myModDiffL]
J = 30.
beta=1.2
gamma=9.5
for model in models:
    if model.__name__ == 'KVORLowerK':
        args = (f0, K)
    elif model.__name__ == 'myModLowerK':
        args = (K,)
    elif model.__name__ == 'myModDiffL':
        args = (beta, gamma)
    else:
        args=()
            
    C = model(*args)
    C.set_xs(np.array([0., 0., -28., 30., 30., 30., -15., -15.]))
    wr = Wrapper(C)
    wr.solve(f0=C.f0, J0=J)
    
#     wr.testHyperBind()
    
    folderName = join('/home/const/Dropbox/Documents/For DN/Very final/data',
                      model.__name__+'J=%.1f'%(J))
    if model.__name__ == 'KVORLowerK':
        folderName += '%.2f_%.0f'%(f0, K)
        
    if model.__name__ == 'myModLowerK':
        folderName += '%.0f'%(K)        
#     wr.dumpEos(folderName)
#     exit()
    print folderName
    
#     wr.dumpMasses(folderName)
#     wr.dumpEos(folderName)
    wr.dumpAll(folderName, folderName+'.zip')
#     wr.dumpScalings(folderName)
#     wr.dumpChi(folderName)
#     wr.dumpMassesCrust(ncut_crust=0.6, ncut_eos=0.8, folderName=folderName)
#     wr.dumpLandauParams(folderName)


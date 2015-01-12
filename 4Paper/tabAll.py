import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join


# models = [Models.KVOR, Models.KVOR_cut_mod, Models.myMod]
models = [Models.myMod]
# models=[Models.KVOR_cut_02]
# models = [Models.KVOR_cut_05]
K = 240
f0 = 0.27
# models = [Models.KVOR_cut_0196_narrow]
# models = [Models.KVOR_cut_02, Models.KVOR_cut_03, Models.KVOR_cut_mod,
#           Models.KVOR_cut_05]

for model in models:
    if model.__name__ == 'KVORLowerK':
        args = (f0, K)
    elif model.__name__ == 'myModLowerK':
        args = (K,)
    else:
        args=()
            
    C = model(*args)
    C.set_xs(np.array([0., 0., -28., 30., 30., 30., -15., -15.]))
    wr = Wrapper(C)


#     wr.testHyperBind()
    
    folderName = join('/home/const/Dropbox/Documents/For DN/Very final/data',
                      model.__name__)
    if model.__name__ == 'KVORLowerK':
        folderName += '%.2f_%.0f'%(f0, K)
        
    if model.__name__ == 'myModLowerK':
        folderName += '%.0f'%(K)        
#     wr.dumpEos(folderName)
#     exit()
    print folderName
    
#     C.Hyper=0
#     wr.dumpMasses(folderName)
#     wr.dumpEos(folderName)
#     wr.dumpAll(folderName, folderName+'.zip')
#     wr.dumpHyper(folderName, npoints=400, verbose=1)
#     wr.dumpScalings(folderName)
#     wr.dumpChi(folderName)
#     wr.dumpUofE(folderName)
#     wr.dumpUofE_anti(folderName)
#     wr.dumpLandauParams(folderName)
#     wr.dumpLandauParamsNM(folderName)
#     wr.dumpVs(folderName)
    C.Hyper = 0
    C.phi_meson=1
# 
#     wr.dumpJ(folderName)
    wr.dumpMassesCrust(ncut_crust=0.6, ncut_eos=0.8, folderName=folderName)
#     wr.dumpLandauParams(folderName)


import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join


# models = [Models.KVOR, Models.KVOR_cut_mod, Models.myMod]
# models = [Models.myMod3]
models = [Models.KVOR_cut_mod]
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
    
#     wr.dumpMasses(folderName)
#     wr.dumpEos(folderName)
#     wr.dumpAll(folderName, folderName+'.zip')
    wr.dumpScalings(folderName)
#     wr.dumpMassesCrust(ncut_crust=0.6, ncut_eos=0.8, folderName=folderName)
#     wr.dumpLandauParams(folderName)


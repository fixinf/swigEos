import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join


models = [Models.waleckaMatsuiCut]

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
    wr.dumpAll(folderName, folderName+'.zip')
    wr.dumpScalings(folderName)
    wr.dumpChi(folderName)
    wr.dumpUofE(folderName)
    wr.dumpUofE_anti(folderName)
    wr.dumpLandauParams(folderName)
    wr.dumpLandauParamsNM(folderName)
    wr.dumpVs(folderName)
    C.Hyper = 0
    C.phi_meson=1

#     wr.dumpJ(folderName)
    wr.dumpMassesCrust(ncut_crust=0.6, ncut_eos=0.8, folderName=folderName)
#     wr.dumpLandauParams(folderName)

import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join


models = [Models.KVOR, Models.KVOR_cut_mod, Models.KVOR_cut_02, Models.KVOR_cut_03, Models.myMod]
# models = [Models.myModNoPhi]
# models=[Models.KVOR_cut_mod]
# models = [Models.KVOR_cut_05]
K = 240
f0 = 0.27
# models = [Models.KVOR_cut_0196_narrow]
# models = [Models.myModNoPhi]#, Models.KVOR_cut_mod]
#           Models.KVOR_cut_05]

# models = [Models.KVOR_cut_mod, Models.KVOR_cut_03, Models.KVOR_cut_02]

# for model in models:
#     m = model()
#     print model.__name__
#     print m.Cs, m.Co, m.Cr, m.b, m.c
#     print m.omega_a, m.omega_b, m.omega_f

DUs = []
# exit()
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
    C.Hyper=0
#     wr.showVsSymLow()
#     wr.dumpMasses(folderName)
#     wr.dumpEos(folderName)
#     wr.dumpAll(folderName, folderName+'.zip')
#     wr.dumpScalings(folderName)
#     wr.dumpChi(folderName)
#     wr.dumpUofE(folderName)
#     wr.dumpUofE_anti(folderName)
#     wr.dumpLandauParams(folderName)
#     wr.dumpLandauParamsNM(folderName)
#     wr.dumpVs(folderName)
    C.Hyper = 0
    nrange = np.linspace(0., 0.005, 100)
    f0 = wr.f0(nrange)
    plt.plot(nrange/wr.n0, 1+f0)
    plt.show()  
#   
#     wr.dumpJ(folderName)
#     n,m,r,mg1,mg2 = wr.dumpMassesCrust(ncut_crust=0.6, ncut_eos=0.8, folderName=folderName)
#     n_max = n[np.argmax(m)]
#     print n_max, max(m)
#     C.Hyper = 0
#     wr.dumpStarProfile(n_max, folderName=folderName)
#     wr.dumpLandauParams(folderName)
#     DUs.append(wr.eval_DU(crust=1))
#     C.sigma_kind = 0
#     wr.dumpHyper(folderName, npoints=400, verbose=1)
#     
#     if model.__name__ == 'myMod':
#         C.sigma_kind = 1
#         print '!!!!!!!!!!!!!!!!!!!!!!!!!!'
#         zl = 3.
#         zx = 3.
#         alpha = 2.5
#         C.phi_meson = 0
#         C.set_hs_alpha(np.array([0., 0.] + [alpha for i in xrange(6)]))
#         C.set_hs_z(np.array([0., 0., zl, 0., 0., 0., zx, zx]))
#         wr.dumpHyper(folderName+'Sigma', npoints=400, verbose=1)
#         
#         zl = 1.
#         zx = 1.
#         alpha = 1.
#         C.phi_gamma = 1.
#         C.phi_z = 1
#         C.phi_meson = 1
#         
#         C.set_hs_alpha(np.array([0., 0.] + [alpha for i in xrange(6)]))
#         C.set_hs_z(np.array([0., 0., zl, 0., 0., 0., zx, zx]))
#         wr.dumpHyper(folderName+'SigmaPhi', npoints=400, verbose=1)
    
print DUs
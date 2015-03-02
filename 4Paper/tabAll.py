import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join


models = [Models.KVOR, Models.KVOR_cut_mod, Models.KVOR_cut_02, Models.KVOR_cut_03, Models.myMod]
# models = [Models.KVOR, Models.KVOR_cut_03, Models.myMod, Models.KVOR_cut_03Chi, Models.myModChi]
# models = [Models.myModNoPhi]
# models = [Models.myModChi]
# models=[Models.myModChi]
# models = [Models.KVOR, Models.KVOR_cut_02]
K = 240
f0 = 0.27
# models = [Models.KVOR_cut_03_nosolve]
# models = [Models.KVOR_cut_03]
# models = [Models.myModNoPhi]#, Models.KVOR_cut_mod]
#           Models.KVOR_cut_05]
# for model in models:
#     C = model()
#     print C.Co/C.eta_o(C.f0)
# 
# exit()
# models = [Models.KVOR_cut_mod, Models.KVOR_cut_03, Models.KVOR_cut_02]
# models = [Models.waleckaMatsui, Models.waleckaMatsuiCut]
# for model in models:
#     m = model()
#     print model.__name__
#     print m.Cs, m.Co, m.Cr, m.b, m.c
#     print m.omega_a, m.omega_b, m.omega_f

DUs = []
# exit()
# for model in models:
#     C = model()
#     wr = Wrapper(C)
#     wr.testHyperBind()
# exit()
xss = []
for model in models:
    C = model()
    xss.append([model.__name__, [C.X_s[i] for i in range(8)]])
    
for i in xss:
    print i
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
#     exit()

#     n = np.linspace(0, 32*wr.n0, 200)
#     E, f = wr.Esymm(n, ret_f=1)
#     plt.plot(n/wr.n0, f)
#     plt.show()
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
# 
    wr.dumpJ(folderName)
# #     exit()
#     C.Hyper = 0
#     wr.dumpVs(folderName)
#     wr.dumpVsSymLow(folderName)
#    
#     wr.dumpJ(folderName)
#     n,m,r,mg1,mg2 = wr.dumpMassesCrust(ncut_crust=0.45, ncut_eos=0.6, folderName=folderName, show=0, inter='cubic')
#     n_max = n[np.argmax(m)]
#     print n_max, max(m)
#     C.Hyper = 0
#     wr.dumpStarProfile(n_max, folderName=folderName)
#     wr.dumpLandauParams(folderName)
#     DUs.append(wr.eval_DU(crust=1))
#     C.sigma_kind = 0

#     wr.dumpHyper(folderName, npoints=1000, verbose=1)
#     wr.dumpVsHyper(folderName)
#     wr.dumpHyperScalings(folderName)

#     if ((model.__name__=='KVOR_cut_03Chi')
#         or (model.__name__ == 'myModChi')
#         ):
#         Hyper = 1
#         C.sigma_kind = 1
# #         print '!!!!!!!!!!!!!!!!!!!!!!!!!!'
# #         zl = 3.
# #         zx = 3.
# #         alpha = 2.5
# #         C.phi_meson = 0
# #         C.set_hs_alpha(np.array([0., 0.] + [alpha for i in xrange(6)]))
# #         C.set_hs_z(np.array([0., 0., zl, zl, zl, zl, zx, zx]))
# # #         C.set_hs_z(np.array([0., 0., zl, 0., 0., 0., zx, zx]))
# #         wr.dumpHyper(folderName+'Sigma', npoints=400, verbose=1)
# #         wr.dumpHyperScalings(folderName+'Sigma')
#         C.phi_a = 0.
#         C.hyper_sigma_kind = 0
#         app='0'
#          
#         zl = 2.
#         zx = 2.
#         alpha = 2.
#         C.phi_gamma = 1.
#         C.phi_z = 1
#         C.phi_meson = 1
#         C.phi_kind = 1 
#         C.set_hs_alpha(np.array([0., 0.] + [alpha for i in xrange(6)]))
# #         C.set_hs_z(np.array([0., 0., zl, 0., 0., 0., zx, zx]))
#         C.set_hs_z(np.array([0., 0., zl, zl, zl, zl, zx, zx]))
#         wr.dumpHyper(folderName+'SigmaPhi'+app, npoints=400, verbose=1)
# #         wr.dumpHyperScalings(folderName+'SigmaPhi'+app)
# #         wr.dumpVsHyper(folderName+'SigmaPhi'+app)
    
print DUs
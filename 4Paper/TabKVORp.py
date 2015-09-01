import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join



# models=[Models.myMod]
models = [Models.KVORphi]

DUs = []

xss = []
for model in models:
    C = model()
    xss.append([model.__name__, [C.X_s[i] for i in range(8)]])
    
for i in xss:
    print(i)
# exit()

for model in models:

            
    C = model()
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
#     print C.eta_r(0.6146)
#     exit()
#     wr.dumpEos(folderName)
#     exit()
#     wr.testDanielewicz()
#     wr.dumpVs()
#     n = np.linspace(0, 8*wr.n0, 100)
#     f0 = wr.f0(n)
#     plt.plot(n/wr.n0, f0)
#     plt.show()
    print(folderName)
    C.Hyper=0
#     wr.showVsSymLow()
#     wr.dumpMasses(folderName)
#     wr.dumpEos(folderName)
 
 
#     wr.dumpAll(folderName, folderName+'.zip')
#     wr.dumpScalings(folderName)
#     wr.dumpChi(folderName)
#     wr.dumpUofE(folderName)
#     wr.dumpUofE_anti(folderName)
#     C.phi_meson = 0
#     wr.dumpUn(folderName)
#     wr.dumpLandauParams(folderName)
#     wr.dumpScalingsSym(folderName)
#          
#     wr.dumpMassesCrustHyper(folderName=folderName, ret_str=0)
#     C.phi_meson = 0
#     wr.dumpLandauParams(folderName)
#     wr.dumpLandauParamsNM(folderName)
# # 
#     wr.dumpJ(folderName)
# # #     exit()
# #     C.Hyper = 0
#     wr.dumpVs(folderName)
#     wr.dumpVsSymLow(folderName)
# #    
# #     wr.dumpJ(folderName)
#     n,m,r,mg1,mg2 = wr.dumpMassesCrust(ncut_crust=0.45, ncut_eos=0.6, folderName=folderName, show=0, inter='cubic')
#     n_max = n[np.argmax(m)]
#     print n_max, max(m)
# #     C.Hyper = 0
#     wr.dumpStarProfile(n_max, folderName=folderName)
# #     wr.dumpLandauParams(folderName)
#     DUs.append(wr.eval_DU(crust=1))
# # #     C.sigma_kind = 0
# #     
    wr.dumpHyper(folderName, npoints=1000, verbose=1)
#     wr.dumpVsHyper(folderName)
#     wr.dumpHyperScalings(folderName)
#     wr.dumpMassesCrustHyper(folderName=folderName, ret_str=1)


    if ((model.__name__=='KVOR_cut_03Chi')
        or (model.__name__ == 'myModChi')
        or (model.__name__ == 'KVORphi')
        ):
        Hyper = 1
        C.sigma_kind = 1
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!')
        zl = 3.
        zx = 3.
        alpha = 2.5
        C.phi_meson = 0
        C.set_hs_alpha(np.array([0., 0.] + [alpha for i in range(6)]))
        C.set_hs_z(np.array([0., 0., zl, zl, zl, zl, zx, zx]))
#         C.set_hs_z(np.array([0., 0., zl, 0., 0., 0., zx, zx]))
#         wr.dumpHyper(folderName+'Sigma', npoints=400, verbose=1)
#         wr.dumpHyperScalings(folderName+'Sigma')
#         wr.dumpMassesCrustHyper(folderName=folderName+'Sigma', ret_str=1)
        
        C.phi_a = 0.
        C.hyper_sigma_kind = 0
        app='0'
          
        zl = 2.
        zx = 2.
        alpha = 2.
        C.phi_gamma = 1.
        C.phi_z = 1
        C.phi_meson = 1
        C.phi_kind = 1 
        C.set_hs_alpha(np.array([0., 0.] + [alpha for i in range(6)]))
#         C.set_hs_z(np.array([0., 0., zl, 0., 0., 0., zx, zx]))
        C.set_hs_z(np.array([0., 0., zl, zl, zl, zl, zx, zx]))
#         wr.dumpHyper(folderName+'SigmaPhi'+app, npoints=400, verbose=1)
# #         wr.dumpMassesCrustHyper(folderName=folderName+'SigmaPhi'+app, ret_str=1)
#         wr.dumpHyperScalings(folderName+'SigmaPhi'+app)
#         wr.dumpVsHyper(folderName+'SigmaPhi'+app)
#     
print(DUs)
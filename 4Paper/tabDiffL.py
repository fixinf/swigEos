import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join
from time import sleep
from scipy.misc.common import derivative
from scipy import optimize


# models = [Models.KVOR, Models.KVOR_cut_mod, Models.myMod]
# models = [Models.myMod3]
models = [Models.myMod240_L]
gamma = 9
beta = 1.4
J=32.

plist = [[0.8, 7.5],
         [1.2, 7.5],
         [1.5, 7.5],
         [1.4, 9.]
         ]

def f(x, C, C1):
    C.rho_a, C.rho_sat_val,C.beta, C.rho_a1, C.rho_a0, C.rho_d = x
    rho_f = C.rho_f
    res =[C.eta_r(0.27)-1, derivative(C.eta_r, 0.27, dx=1e-3) - 
            derivative(C1.eta_r, 0.27, dx=1e-3), C.eta_r(C.rho_f) - C1.eta_r(rho_f)]
    res.append(C1.eta_r(0) - C.eta_r(0))
    res.append(C1.eta_r(0.1) - C.eta_r(0.1))
#     res.append(C1.eta_r(0.5) - C.eta_r(0.5))
    res.append(C1.eta_r(0.6) - C.eta_r(0.6))
    res = np.array(res)
    return np.sum(res*res)
    

Jlist = [28., 30., 32.]
Jlist = [30.]
# plist = [[1.4, 9.]]
for J in Jlist:
    for p in plist:
        for model in models:
            beta = p[0]
            gamma = p[1]
            C = model(beta=beta, gamma=gamma)
            C.Hyper = 0
#             C.Cr = 70.2979903917
            C.set_xs(np.array([0., 0., -28., 30., 30., 30., -15., -15.]))
            wr = Wrapper(C)
            for j in np.linspace(wr.J(), J, 5):
                wr.solve(f0=C.f0, J0=j, K0=wr.K())
            print wr.L()
    #         exit()
        #     wr.testHyperBind()
            
#             C = Models.myMod()
#             C.rho_f = 0.75
# #             C.rhO_d = 0
#             res = optimize.minimize(lambda z: f(z, C, C2), [C.rho_a, C.rho_sat_val, C.beta, C.rho_a1, C.rho_a0, C.rho_d])
#             print res
            frange = np.linspace(0, 1, 100)
#             plt.plot(frange, map(C.eta_r, frange), frange, map(C2.eta_r,frange))
#             plt.show()
#             wr = Wrapper(C)
#             wr.reset(npoints=100)
#             wr2.reset(npoints=100)
            wr.reset(npoints=100)
            plt.plot(wr.n/wr.n0, wr.concentrations())
            plt.plot(wr.n/wr.n0, wr.rho[:,0])
            plt.show()
            plt.plot(frange, map(C.eta_r, frange))
            plt.show()
            folderName = join('/home/const/Dropbox/Documents/For DN/Very final/data',
                              model.__name__+
                              'J=%.0f_L=%.2f'%(J, wr.L()))
        #     wr.dumpEos(folderName)
        #     exit()
            print folderName
            
        #     wr.dumpMasses(folderName)
#             wr.dumpEos(folderName)
        #     wr.dumpAll(folderName, folderName+'.zip')
        #     wr.dumpScalings(folderName)
#             wr.dumpMassesCrust(ncut_crust=0.6, ncut_eos=0.8, folderName=folderName)
        #     wr.dumpLandauParams(folderName)
            print wr.L()
            wr.dumpLandauParamsNM(folderName)
#             wr.dumpJ(folderName)

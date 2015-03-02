import eosWrap as eos
import matplotlib
from math import pi
matplotlib.use('QT4Agg')
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join
from scipy.misc.common import derivative
from scipy import optimize

KV = Models.KVOR()

models = [Models.waleckaMatsuiCut]
C0 = Models.waleckaMatsui()
wr0 = Wrapper(C0)

# print C0.n0
# print C0.n0*(135./197.33)**3
# print 135/(197.33) * (3*pi**2*KV.n0/2)**(1./3.)
# print 135/(197.33) * (3*pi**2*C0.n0/2)**(1./3.)
# print ((1.42*197.33/135)**3) / ((3*pi**2)/2.)
# print wr0.n0
# print wr0.J(), wr0.K(), wr0.L(), wr0.Ksymm(), wr0.Kprime()
# n = [0.2, wr0.n0]
# es, f = wr0.Esymm( n, ret_f=1)
# print es, f
# print wr0.ESbind(n)
# print eos.f_eq(np.array([C0.n0/2, C0.n0/2]), np.array([0.]), 1, C0)
# 
# print wr0.eval_DU()
# # print wr.eval_DU()
# 
# exit()


C = models[0]()
frange = np.linspace(0., 1, 100) 
wr = Wrapper(C)

C2 = Models.KVOR_cut_mod()
wr2 = Wrapper(C2)
plt.plot(frange, map(C.eta_o, frange), frange, map(C2.eta_o, frange))
plt.show()

C2.alpha = 0
C2.f0 = C.f0

C2.omega_f = 0.55

plt.plot(frange, map(C.eta_o, frange), frange, map(C2.eta_o, frange))
plt.show()
n = np.linspace(0., 2., 100)
C2.beta1 = 0
C2.beta2 = 0
C2.d2om = 0
C2.dom = 0
C2.drho = 0
C2.rho = 0
C2.d2rho = 0
C2.gamma2= 0
C2.rho_a = 0
C2.rho_a2 = 0
C2.rho_a = 0
C2.c1 = 0
C2.rho_a0 = 0
C2.rho_a1 = 0
C2.rho_a2 = 0
C2.rho_a3 = 0
C2.rho_a4 = 0
C2.rho_e = 0
C2.rho_b_low = 0
C2.rho_sat_a = 0
C2.rho_sat_f1 = 1
C2.rho_sat_f2 = 1
C2.rho_tan_a = 0
C2.rho_tan_b = 0
C2.rho_tan_c = 0
C2.omega_b = 0
C2.rho_sat_val = 0
C2.rho_width_f = 1
C2.rho_width_power = 0
C2.gamma = 0
C2.beta = 0
C2.rho_a = 0
def f(x):
    C2.omega_b = x[0]
    C2.omega_f = x[1]
    fm = 0.54
    res = np.array([C.eta_o(fm) - C2.eta_o(fm),
                    derivative(C.eta_o, fm, dx=1e-3) - 
                    derivative(C2.eta_o, fm, dx=1e-3)])
    return np.sum(res**2)
    
res = optimize.minimize(f, [C.omega_b, C.omega_f])
print res
if res.x[0] > 10**5:
    pass

# x_Matsui = [53.30066647,   0.56083489]

plt.plot(frange, map(C.eta_o, frange), frange, map(C2.eta_o, frange))
plt.show()

wr2.n0 = wr.n0

C2.rho_kind = 9
C2.b = 0
C2.c = 0
C2.Cs = C.Cs
C2.Co = C.Co
C2.Cr = C.Cr
plt.plot(n, wr.ESbind(n), n, wr2.ESbind(n))
plt.show()
plt.plot(frange, map(C2.eta_s, frange))
plt.plot(frange, map(C2.eta_r, frange))
plt.plot(frange, map(C2.U, frange))
plt.ylim([0, 2])
plt.show()

C2.d = 0
es0, fs0 = wr.Eneutr(n, ret_f=1)
es2, fs2 = wr2.Eneutr(n, ret_f=1)

plt.plot(n, fs0, n, fs2)
plt.show()



wr0.reset(npoints=100)
wr.reset(npoints=100)
plt.plot(wr.n/wr.n0, wr0.concentrations(), wr.n/wr.n0, wr.concentrations())
plt.show()

C0.Hyper = 0
C.Hyper = 0
# wr0.eval_DU(crust=1)
# wr.eval_DU(crust=1)
# wr0.testDanielewicz()
# wr.testDanielewicz()
# wr0.testPodsiedlowski()
# wr.testPodsiedlowski()
models = [Models.waleckaMatsui]
for model in models:

    C = model()
    C.Hyper=0
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
#     wr.dumpScalings(folderName)
#     wr.dumpChi(folderName)
#     wr.dumpUofE(folderName)
#     wr.dumpUofE_anti(folderName)
#     wr.dumpLandauParams(folderName) 
#     wr.dumpLandauParamsNM(folderName)
#     wr.dumpVs(folderName)
    C.Hyper = 0
    C.phi_meson=1

#     wr.dumpJ(folderName)
    wr.dumpMassesCrust(ncut_crust=0.6, ncut_eos=0.8, folderName=folderName)
#     wr.dumpLandauParams(folderName)

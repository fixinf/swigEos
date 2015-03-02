import matplotlib.pyplot as plt
from Wrapper import Wrapper
import eosWrap as eos
import numpy as np
from scipy.misc.common import derivative
from scipy import optimize
from math import pi



def KVOR():
    C = eos.KVOR()
    C.Csp = 1.
    C.Cs = 179.56233875171566
    C.Co =87.5996397368707
    C.Cr = 100.63642242798424
    C.b = 0.00773460805148428
    C.c = 0.00034461786646922604
    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def waleckaMatsui():
    C= eos.Walecka()
    C.Csp = 1.
    C.SetHyperConstants(2)
    C.Cs = 266.9
    C.Co = 195.7
    C.Cr = 54.71
    C.b = 0
    C.c = 0
    C.n0 = ((1.42*197.33/135)**3) / ((3*pi**2)/2.)
    C.f0 = 0.44
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def KVOR_cut():
    C = eos.KVOR_mod2()
    C.Csp = 1.
    C.Cs = 179.5623289954
    C.Co =87.5996301756
    C.Cr = 100.6364192718 
    C.b = 0.0077346088
    C.c = 0.0003446263
    C.SetHyperConstants(2)
    C.omega_kind = 1
    C.omega_a = 1000
    C.omega_f = 0.4
    C.phi_gamma = 3.
    C.phi_z = 3.5

    wr = Wrapper(C)
    wr.solve(f0=C.f0)
    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def KVOR_cut_mod():
    C = eos.KVOR_cut()
    C.SetHyperConstants(2)
    C.omega_a = 300.
    C.omega_f = 0.3
    C.omega_kind = 1
    C.Csp = 1.
#   
    C.Cs = 179.5623289954
    C.Co = 87.5996301756
    C.Cr = 100.6364192718 
    C.b = 0.0077346088
    C.c = 0.0003446263
    C.omega_a = 300.
    C.omega_f = 0.4
#     C.rho_f = 0.3
#     C.rho_a = 1000.
    
    
    C.phi_gamma = 3.
    C.phi_z = 3.5
    wr = Wrapper(C)
    wr.solve(f0=C.f0)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
#     return C
    return KVOR_tan_04(C)

def myModChi():
    C = myMod()
    C.phi_kind = 1
    return C

def f_cut(x, C, C1):
    C.omega_b = x[0]
    C.omega_f = x[1]
    res = [C.eta_o(C.omega_f) - C1.eta_o(C.omega_f),
            derivative(C.eta_o, C.omega_f, dx=1e-3) - derivative(C1.eta_o, C.omega_f, dx=1e-3)]
    res = np.array(res)
    return np.sum(res*res)

def KVOR_tan_04(C1):
    C = eos.KVOR_cut()
    C.SetHyperConstants(2)
    C.omega_a = 300.
    C.omega_f = 0.3
    C.omega_kind = 1
    C.Csp = 1.
#   
    C.Cs = 179.5623289954
    C.Co = 87.5996301756
    C.Cr = 100.6364192718 
    C.b = 0.0077346088
    C.c = 0.0003446263
    C.omega_a = 300.
    C.omega_f = 0.4
#     C.rho_f = 0.3
#     C.rho_a = 1000.
    
    C.omega_kind = 2
    C.omega_b = 100
    C.omega_a = -0.5
    C.omega_f += 0.05
    
#     C1 = KVOR_cut_mod()
    
    res = optimize.minimize(lambda z: f_cut(z, C, C1), [C.omega_b, C.omega_f])
    
    C.phi_gamma = 3.
    C.phi_z = 3.5
    wr = Wrapper(C)
    wr.solve(f0=C.f0)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def KVOR_cut_03():
#     return KVOR_tan_03()
    C = eos.KVOR_cut()
    C.SetHyperConstants(2)
    C.omega_kind = 1
    C.Csp = 1.
    C.Cs = 179.5623289954
    C.Co =87.5996301756
    C.Cr = 100.6364192718 
    C.b = 0.0077346088
    C.c = 0.0003446263

# 
    C.omega_a = 200.
    C.omega_f = 0.3
    C.rho_f = 0.2
    C.rho_a = 200.

    C.phi_gamma = 3.
    C.phi_z = 3.5
    wr = Wrapper(C)
    wr.solve(f0=C.f0)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
#     return C
    return KVOR_tan_03(C)

def KVOR_cut_03_nosolve():
#     return KVOR_tan_03()
    C = eos.KVOR_cut()
    C.SetHyperConstants(2)
    C.omega_kind = 1
    C.Csp = 1.
    C.Cs = 179.5623289954
    C.Co =87.5996301756
    C.Cr = 100.6364192718 
    C.b = 0.0077346088
    C.c = 0.0003446263

# 
    C.omega_a = 200.
    C.omega_f = 0.3
    C.rho_f = 0.2
    C.rho_a = 200.

    C.phi_gamma = 3.
    C.phi_z = 3.5
    wr = Wrapper(C)
#     wr.solve(f0=C.f0)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
#     return C
    return KVOR_tan_03(C,solve=0)

def KVOR_cut_03Chi():
    C = KVOR_cut_03()
    C.phi_kind = 1
    C.phi_a = 0
    return C

def KVOR_tan_03(C1, solve=1):
    C = eos.KVOR_cut()
    C.SetHyperConstants(2)
    C.omega_a = 300.
    C.omega_f = 0.3
    C.omega_kind = 1
    C.Csp = 1.
#   
    C.Cs = 179.5623289954
    C.Co = 87.5996301756
    C.Cr = 100.6364192718 
    C.b = 0.0077346088
    C.c = 0.0003446263
    C.omega_a = 300.
    C.omega_f = 0.4
#     C.rho_f = 0.3
#     C.rho_a = 1000.
    
    C.omega_kind = 2
    C.omega_b = 100
    C.omega_a = -0.5
    C.omega_f += 0.05
    
#     C1 = KVOR_cut_03()
    
    res = optimize.minimize(lambda z: f_cut(z, C, C1), [C.omega_b, C.omega_f])
    
    print res.x
    
    C.phi_gamma = 3.
    C.phi_z = 3.5
    wr = Wrapper(C)
    if solve:
        wr.solve(f0=C.f0)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def KVR_cut():
    C = eos.KVOR_mod2()

    C.SetHyperConstants(2)
    wr = Wrapper(C)
    
    C.Csp = 1
    
    C.Cs = 179.56233875157545
    C.Co =  87.59963973682763
    C.Cr = 100.63642242792484
    C.b = 0.007734608051455927
    C.c = 0.0003446178665624873
    
    C.rho_f = 0.48
    C.rho_a = 0*1000.
    C.omega_f = 0.3
    C.omega_a = 1000.
    C.omega_kind = 1
    C.rho_power = 1.
    f0 = 0.2
    C.rho_kind = 1
    
    C.beta = 2.9
    C.gamma = 2.
    C.alpha = 0.0
    C.phi_gamma = 3.
    C.phi_z = 3.5
    C.SetHyperConstants(2)
    
    C.f0 = f0
    wr.solve(f0 = C.f0, E0=-15.8, K0 = 250, J0=28.0, iter = 3000)
    return C

def KVOR_cut_smooth():
    C = eos.KVOR_mod2()
    C.Csp = 1.
    C.Cs = 178.9560101100
    C.Co =87.5996399401
    C.Cr = 100.6364207940
    C.b = 0.0073248882
    C.c = 0.0028050590
    print C.eta_o(0.2), C.eta_p(0.2), C.eta_r(0.2), C.eta_s(0.2)
    C.d = 0.
    C.SetHyperConstants(2)
    C.omega_kind = 1
    C.omega_a = 200
    C.omega_f = 0.387
    C.rho_kind = 1
    C.rho_power = 2
    C.gamma = 5.2
    C.beta = 1.2
    C.alpha = 1.
    
    C.phi_gamma = 3.
    C.phi_z = 3.5
    wr = Wrapper(C)
    wr.solve(f0=C.f0)
    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def KVOR_cut_02():
    C = eos.KVOR_cut()
    C.SetHyperConstants(2)
    C.omega_kind = 1
    C.Csp = 1.
    C.Cs = 179.5623289954
    C.Co =87.5996301756
    C.Cr = 100.6364192718 
    C.b = 0.0077346088
    C.c = 0.0003446263

# 
    C.omega_a = 200.
    C.omega_f = 0.2
    C.rho_f = 0.197
    C.rho_a = 2000.

    C.phi_gamma = 3.
    C.phi_z = 3.5
    wr = Wrapper(C)
    wr.solve(f0=C.f0)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return KVOR_tan_02(C)

def KVOR_cut_02_nosolve():
    C = eos.KVOR_cut()
    C.SetHyperConstants(2)
    C.omega_kind = 1
    C.Csp = 1.
    C.Cs = 179.5623289954
    C.Co =87.5996301756
    C.Cr = 100.6364192718 
    C.b = 0.0077346088
    C.c = 0.0003446263

# 
    C.omega_a = 200.
    C.omega_f = 0.2
    C.rho_f = 0.197
    C.rho_a = 2000.

    C.phi_gamma = 3.
    C.phi_z = 3.5
    wr = Wrapper(C)
#     wr.solve(f0=C.f0)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return KVOR_tan_02(C, solve=False)

def KVOR_tan_02(C1, solve=True):
    C = eos.KVOR_cut()
    C.SetHyperConstants(2)
    C.omega_a = 300.
    C.omega_f = 0.3
    C.omega_kind = 1
    C.Csp = 1.
#   
    C.Cs = 184.2559948015
    C.Co = 87.5937257485
    C.Cr = 100.6364192718 
    C.b = 0.0099933721
    C.c = -0.0076379652
    C.omega_a = 300.
    C.omega_f = 0.4
#     C.rho_f = 0.3
#     C.rho_a = 1000.
    
    C.omega_kind = 2
    C.omega_b = 100
    C.omega_a = -0.2
    C.omega_f = 0.28
    
    res = optimize.minimize(lambda z: f_cut(z, C, C1), [C.omega_b, C.omega_f])
    
    print res.x
    
    C.phi_gamma = 3.
    C.phi_z = 3.5
    wr = Wrapper(C)
#     for f in np.linspace(0.3, C.omega_f, 10):
#         C.omega_f = f
    if solve:
        wr.solve(f0=C.f0)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def KVOR_cut_0196():
    C = eos.KVOR_mod2()
    C.Csp = 1.
    C.Cs = 178.9560101100
    C.Co =87.5996399401
    C.Cr = 100.6364207940
    C.b = 0.0073248882
    C.c = 0.0028050590
    print C.eta_o(0.2), C.eta_p(0.2), C.eta_r(0.2), C.eta_s(0.2)
    C.d = 0.
    C.SetHyperConstants(2)
    C.omega_kind = 1
#     C.omega_a = 1000
#     C.omega_f = 0.4
#Smooth version:
    
    C.rho_kind = 1
    C.rho_power = 2
    C.gamma = 5.2
    C.beta = 1.2
    C.alpha = 1.
    
    C.omega_a = 200
    C.omega_f = 0.196
    
    C.phi_gamma = 3.
    C.phi_z = 3.5
    wr = Wrapper(C)
#     wr.solve(f0=C.f0)
    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def KVOR_cut_0196_narrow():
    C = eos.KVOR_mod2()
    C.Csp = 1.
    C.Cs = 178.9560101100
    C.Co =87.5996399401
    C.Cr = 100.6364207940
    C.b = 0.0073248882
    C.c = 0.0028050590
    print C.eta_o(0.2), C.eta_p(0.2), C.eta_r(0.2), C.eta_s(0.2)
    C.d = 0.
    C.SetHyperConstants(2)
    C.omega_kind = 1
#     C.omega_a = 1000
#     C.omega_f = 0.4
#Smooth version:
    
    C.rho_kind = 1
    C.rho_power = 2
    C.gamma = 5.2
    C.beta = 1.2
    C.alpha = 1.
    
    C.omega_a = 20000
    C.omega_f = 0.196
    
    C.phi_gamma = 3.
    C.phi_z = 3.5
    wr = Wrapper(C)
#     wr.solve(f0=C.f0)
    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return KVOR_tan_0196(C)
#     return C

def KVOR_tan_0196(C1):
    C = eos.KVOR_cut()
    C.SetHyperConstants(2)
    C.omega_a = 300.
    C.omega_f = 0.3
    C.omega_kind = 1
    C.Csp = 1.
#   
    C.Cs = 184.2559948015
    C.Co = 87.5937257485
    C.Cr = 100.6364192718 
    C.b = 0.0099933721
    C.c = -0.0076379652
    C.omega_a = 300.
    C.omega_f = 0.4
#     C.rho_f = 0.3
#     C.rho_a = 1000.
    
    C.omega_kind = 2
    C.omega_b = 120
    C.omega_a = -0.8
    C.omega_f = 0.197
    
    res = optimize.minimize(lambda z: f_cut(z, C, C1), [C.omega_b, C.omega_f])
    
    print res
    
    frange = np.linspace(0, 1, 1000)
    plt.plot(frange, map(C.eta_o, frange))
    plt.plot(frange, map(C1.eta_o, frange))
    plt.show()
    
    C.phi_gamma = 3.
    C.phi_z = 3.5
    wr = Wrapper(C)
#     for f in np.linspace(0.3, C.omega_f, 10):
#         C.omega_f = f
    wr.solve(f0=C.f0)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def KVOR_cut_05():
    C = eos.KVOR_cut()
    C.SetHyperConstants(2)
    C.omega_kind = 1
    C.Csp = 1.
    C.Cs = 179.5623289954
    C.Co =87.5996301756
    C.Cr = 100.6364192718 
    C.b = 0.0077346088
    C.c = 0.0003446263

# 
    C.omega_a = 300.
    C.omega_f = 0.5
    C.rho_f = 0.2
    C.rho_a = 0.

    C.phi_gamma = 3.
    C.phi_z = 3.5
    wr = Wrapper(C)
    wr.solve(f0=C.f0)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C


def myMod2():
    C = eos.KVOR_mod2()
    
    C.Cs = 227.8259926550
    C.Co = 134.8826104284
    C.Cr = 93.1990895430
    C.b = 0.0055916190
    C.c = -0.0121228812 
    C.alpha = 0.85
    C.z = 0.65
    
    C.omega_kind = 1
    
    C.omega_a = 200
    C.omega_f = 0.53
    
    C.beta = 0.8
    C.gamma = 7.5
              
    C.phi_a = -0.85
    C.phi_f = 0.28
    
    C.d = -5.5
    
    C.rho_f = 0.75
    C.rho_a = 1000

    C.f0 = 0.27
    C.rho_kind = 1
    C.rho_power = 2.0
    C.omega_c = -15000
    
    C.SetHyperConstants(2)
    
    C.Csp = 380.
    
    C.phi_gamma = 3.0
    C.phi_z=3.5
    
    
    
    wr = Wrapper(C)
#     wr.solve(f0 = C.f0)
    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def myMod3():
    C = eos.KVOR_mod2()
    
    C.Cs = 227.8259926550
    C.Co = 134.8826104284
    C.Cr = 93.1990895430
    C.b = 0.0055916190
    C.c = -0.0121228812 
    C.alpha = 0.85
    C.z = 0.65
    
    C.omega_kind = 1
    
    C.omega_a = 0
    C.omega_f = 0.695
    
    C.beta = 0.8
    C.gamma = 7.5
              
    C.phi_a = -0.85
    C.phi_f = 0.28
    
    C.d = -5.5
    
    C.rho_f = 0.75
    C.rho_a = 1000

    C.f0 = 0.27
    C.rho_kind = 1
    C.rho_power = 2.0
    C.omega_c = -15000
    
    C.SetHyperConstants(2)
    
    C.Csp = 380.
    
    C.phi_gamma = 3.0
    C.phi_z=3.5
    
    
    
    wr = Wrapper(C)
#     wr.solve(f0 = C.f0)
    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def myModLowerK(K=275.):
    C = eos.KVOR_mod2()
    
    C.Cs = 232.2635684167
    C.Co = 134.8826104284
    C.Cr = 93.1990895430
    C.b = 0.0063104759
    C.c = -0.0119543213  
    C.alpha = 0.85
    C.z = 0.65
    
    C.omega_a = 6.45
    C.omega_f = 0.53
    
    C.beta = 0.8
    C.gamma = 7.5
              
    C.phi_a = -0.85
    C.phi_f = 0.28
    
    C.d = -5.5
    
    C.rho_f = 0.75
    C.rho_a = 1000

    C.f0 = 0.27
    C.rho_kind = 1
    C.rho_power = 2.0
    C.omega_c = -15000
    
    C.SetHyperConstants(2)
    
    C.Csp = 380.
    
    C.phi_gamma = 3.0
    C.phi_z=3.5
    
    wr = Wrapper(C)
    wr.solve(f0 = C.f0, K0=K)
    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def myModDiffL(beta=0.8, gamma = 7.5):
    C = eos.KVOR_mod2()
    
    C.Cs = 227.8259926550
    C.Co = 134.8826104284
    C.Cr = 93.1990895430
    C.b = 0.0055916190
    C.c = -0.0121228812 
    C.alpha = 0.85
    C.z = 0.65
    
    C.omega_a = 6.45
    C.omega_f = 0.53
    
    C.beta = beta
    C.gamma = gamma
              
    C.phi_a = -0.85
    C.phi_f = 0.28
    
    C.d = -5.5
    
    C.rho_f = 0.75
    C.rho_a = 1000

    C.f0 = 0.27
    C.rho_kind = 1
    C.rho_power = 2.0
    C.omega_c = -15000
    
    C.SetHyperConstants(2)
    
    C.Csp = 380.
    
    C.phi_gamma = 3.0
    C.phi_z=3.5
    
    wr = Wrapper(C)
    wr.solve(f0 = C.f0)
    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    C.Hyper = 0
    return C

def myModDiffL_J30(beta=0.8, gamma = 7.5):
    C = eos.KVOR_mod2()
    
    C.Cs = 227.8259926550
    C.Co = 134.8826104284
    C.Cr = 93.1990895430
    C.b = 0.0055916190
    C.c = -0.0121228812 
    C.alpha = 0.85
    C.z = 0.65
    
    C.omega_a = 6.45
    C.omega_f = 0.53
    
    C.beta = beta
    C.gamma = gamma
              
    C.phi_a = -0.85
    C.phi_f = 0.28
    
    C.d = -5.5
    
    C.rho_f = 0.75
    C.rho_a = 1000

    C.f0 = 0.27
    C.rho_kind = 1
    C.rho_power = 2.0
    C.omega_c = -15000
    
    C.SetHyperConstants(2)
    
    C.Csp = 380.
    
    C.phi_gamma = 3.0
    C.phi_z=3.5
    
    wr = Wrapper(C)
    wr.solve(f0 = C.f0, J0=30)
    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def myModDiffL_J28(beta=0.8, gamma = 7.5):
    C = eos.KVOR_mod2()
    
    C.Cs = 227.8259926550
    C.Co = 134.8826104284
    C.Cr = 93.1990895430
    C.b = 0.0055916190
    C.c = -0.0121228812 
    C.alpha = 0.85
    C.z = 0.65
    
    C.omega_a = 6.45
    C.omega_f = 0.53
    
    C.beta = beta
    C.gamma = gamma
              
    C.phi_a = -0.85
    C.phi_f = 0.28
    
    C.d = -5.5
    
    C.rho_f = 0.75
    C.rho_a = 1000

    C.f0 = 0.27
    C.rho_kind = 1
    C.rho_power = 2.0
    C.omega_c = -15000
    
    C.SetHyperConstants(2)
    
    C.Csp = 380.
    
    C.phi_gamma = 3.0
    C.phi_z=3.5
    
    wr = Wrapper(C)
    wr.solve(f0 = C.f0, J0=28.)
    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def KVORLowerK(f0, K):
    C = eos.KVOR_mod2()
    C.Csp = 1.
    C.Cs = 178.9560101100
    C.Co =87.5996399401
    C.Cr = 100.6364207940
    C.b = 0.0073248882
    C.c = 0.0028050590
    print C.eta_o(0.2), C.eta_p(0.2), C.eta_r(0.2), C.eta_s(0.2)
    C.d = 0.
    C.SetHyperConstants(2)
    C.omega_kind = 1
#     C.omega_a = 1000
#     C.omega_f = 0.4
#Smooth version:
    
    C.rho_kind = 1
    C.rho_power = 2
    C.gamma = 5.2
    C.beta = 1.2
    C.alpha = 1.
    
    C.omega_a = 0.
    C.omega_f = 0.387
    
    C.phi_gamma = 3.
    C.phi_z = 3.5
    C.f0 = f0
    wr = Wrapper(C)
    wr.solve(f0=f0, K0=K)
    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def myModLow():
    C = eos.KVOR_mod2()
    
    C.Cs = 227.8259926550
    C.Co = 134.8826104284
    C.Cr = 93.1990895430
    C.b = 0.0055916190
    C.c = -0.0121228812 
    C.alpha = 0.85
    C.z = 0.65
    
    C.omega_a = 6.45
    C.omega_f = 0.53
    
    C.beta = 0.8
    C.gamma = 7.5
              
    C.phi_a = -0.85
    C.phi_f = 0.28
    
    C.d = -5.5
    
    C.rho_f = 0.75
    C.rho_a = 1000

    C.f0 = 0.27
    C.rho_kind = 1
    C.rho_power = 2.0
    C.omega_c = -15000
    
    C.SetHyperConstants(2)
    
    C.Csp = 380.
    
    C.phi_gamma = 3.0
    C.phi_z=3.5
    
    wr = Wrapper(C)
#     wr.solve(f0 = C.f0)
    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C


def myMod240():
    C = eos.KVOR_mod2()
    
    
    C.Cs = 234.1580555799
    C.Co = 134.8826104616
    C.b = 0.0046776700
    C.c = -0.0029781609
    C.Cr = 93.1990895430
    
    C.alpha = 0.4
    C.z = 0.65
    
    C.omega_a = 0.8
    C.omega_f = 0.55
    
    C.beta = 0.8
    C.gamma = 7.5
              
    C.phi_a = -0.0
    C.phi_f = 0.28
    
    C.d = -0.5
    
    C.rho_f = 0.75
    C.rho_a = 1000

    C.f0 = 0.27
    C.rho_kind = 1
    C.rho_power = 2.0
    C.omega_c = -20000
    
    C.SetHyperConstants(2)
    
    C.Csp = 380.
    
    C.phi_gamma = 3.0
    C.phi_z=3.5
    
    wr = Wrapper(C)
    wr.solve(f0 = C.f0, K0=240.)
    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def myMod240_L(beta=0.8, gamma=7.5):
    C = eos.KVOR_mod2()
    
    
    C.Cs = 234.1580555799
    C.Co = 134.8826104616
    C.b = 0.0046776700
    C.c = -0.0029781609
    C.Cr = 81.7485399416
    
    C.alpha = 0.4
    C.z = 0.65
    
    C.omega_a = 0.8
    C.omega_f = 0.55
    
    C.beta = beta
    C.gamma = gamma
              
    C.phi_a = -0.0
    C.phi_f = 0.28
    
    C.d = -0.5
    
    C.rho_f = 0.75
    C.rho_a = 1000

    C.f0 = 0.27
    C.rho_kind = 1
    C.rho_power = 2.0
    C.omega_c = -20000
    
    C.SetHyperConstants(2)
    
    C.Csp = 380.
    
    C.phi_gamma = 3.0
    C.phi_z=3.5
    
    wr = Wrapper(C)
    
    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    C.Hyper = 0
    return C


def myMod240J30():
    C = eos.KVOR_mod2()
    
    
    C.Cs = 234.1580555799
    C.Co = 134.8826104616
    C.b = 0.0046776700
    C.c = -0.0029781609
    C.Cr = 93.1990895430
    
    C.alpha = 0.4
    C.z = 0.65
    
    C.omega_a = 0.8
    C.omega_f = 0.55
    
    C.beta = 0.8
    C.gamma = 7.5
              
    C.phi_a = -0.0
    C.phi_f = 0.28
    
    C.d = -0.5
    
    C.rho_f = 0.75
    C.rho_a = 1000

    C.f0 = 0.27
    C.rho_kind = 1
    C.rho_power = 2.0
    C.omega_c = -20000
    
    C.SetHyperConstants(2)
    
    C.Csp = 380.
    
    C.phi_gamma = 3.0
    C.phi_z=3.5
    
    wr = Wrapper(C)
    wr.solve(f0 = C.f0, K0=240., J0=30.)
    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def myMod240J28():
    C = eos.KVOR_mod2()
    
    
    C.Cs = 234.1580555799
    C.Co = 134.8826104616
    C.b = 0.0046776700
    C.c = -0.0029781609
    C.Cr = 93.1990895430
    
    C.alpha = 0.4
    C.z = 0.65
    
    C.omega_a = 0.8
    C.omega_f = 0.55
    
    C.beta = 0.8
    C.gamma = 7.5
              
    C.phi_a = -0.0
    C.phi_f = 0.28
    
    C.d = -0.5
    
    C.rho_f = 0.75
    C.rho_a = 1000

    C.f0 = 0.27
    C.rho_kind = 1
    C.rho_power = 2.0
    C.omega_c = -20000
    
    C.SetHyperConstants(2)
    
    C.Csp = 380.
    
    C.phi_gamma = 3.0
    C.phi_z=3.5
    
    wr = Wrapper(C)
    wr.solve(f0 = C.f0, K0=240., J0=30.)
    wr.solve(f0 = C.f0, K0=240., J0=29.)
    wr.solve(f0 = C.f0, K0=240., J0=28.)

    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def myMod():
    return myModFinal()
    return myModR()
    C = myMod240_L(1.2, 7.5)
    C.Cs = 234.1580555929
    C.Co = 134.8826104259
    C.Cr = 81.7485393895
    C.b=  0.0046776700
    C.c = -0.0029781609
    wr = Wrapper(C)
#     for j in np.linspace(wr.J(), 30., 4):
#         wr.solve(f0=C.f0, K0=wr.K(), J0=j)
    wr.solve(f0=C.f0, J0=30., K0=240.)
    print wr.L(), wr.J()
#     exit()
    return C
    
    C.Cs = 227.8259926550
    C.Co = 134.8826104284
    C.Cr = 93.1990895430
    C.b = 0.0055916190
    C.c = -0.0121228812 
    C.alpha = 0.85
    C.z = 0.65
    
    C.omega_a = 6.45
    C.omega_f = 0.53
    
    C.beta = 0.8
    C.gamma = 7.5
              
    C.phi_a = -0.85
    C.phi_f = 0.28
    
    C.d = -5.5
    
    C.rho_f = 0.75
    C.rho_a = 1000

    C.f0 = 0.27
    C.rho_kind = 1
    C.rho_power = 2.0
    C.omega_c = -15000
    
    C.SetHyperConstants(2)
    
    C.Csp = 380.
    
    C.phi_gamma = 0.
    C.phi_z = 0.
    
#     C.phi_gamma = 3.0
#     C.phi_z=3.5
    
    wr = Wrapper(C)
#     wr.solve(f0 = C.f0)
    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def myModSat():
    C = myModR()
    C.rho_kind = 3
    C.rho_sat_val = 6.5185
    C.rho_sat_f1 = 0.5;
    C.rho_sat_f2 = 0.;
    C.rho_sat_a = 20
    wr = Wrapper(C)
    wr.solve(f0=C.f0, K0=240., J0=30.)
    return C



def myModR2():
    C = myMod240_L(1.2, 7.5)
    C.Cs = 234.1580555929
    C.Co = 134.8826104259
    C.Cr = 81.7485393895
    C.b=  0.0046776700
    C.c = -0.0029781609
    C.rho_kind = 6
    C.rho_a = 0.6
    C.rho_sat_a = 20
    C.rho_sat_f1 = 0.35
    C.beta = 3.4
    C.gamma = 0
    C.rho_sat_val = 14.58227877
    C.rho_f = 0.65309008
    wr = Wrapper(C)
#     for j in np.linspace(wr.J(), 30., 4):
#         wr.solve(f0=C.f0, K0=wr.K(), J0=j)
    wr.solve(f0=C.f0, J0=30., K0=240.)
    print wr.L(), wr.J()
#     exit()
    return C

def myModR():
    C = myMod240_L(1.2, 7.5)
    C.Cs = 234.1580555929
    C.Co = 134.8826104259
    C.Cr = 81.7485393895
    C.b=  0.0046776700
    C.c = -0.0029781609
    C.rho_f = 0.45
    C.rho_a = 100.
    wr = Wrapper(C)
#     for j in np.linspace(wr.J(), 30., 4):
#         wr.solve(f0=C.f0, K0=wr.K(), J0=j)
    wr.solve(f0=C.f0, J0=30., K0=240.)
    print wr.L(), wr.J()
#     exit()
    return C
    
    C.Cs = 227.8259926550
    C.Co = 134.8826104284
    C.Cr = 93.1990895430
    C.b = 0.0055916190
    C.c = -0.0121228812 
    C.alpha = 0.85
    C.z = 0.65
    
    C.omega_a = 6.45
    C.omega_f = 0.53
    
    C.beta = 0.8
    C.gamma = 7.5
              
    C.phi_a = -0.85
    C.phi_f = 0.28
    
    C.d = -5.5
    
    C.rho_f = 0.75
    C.rho_a = 1000

    C.f0 = 0.27
    C.rho_kind = 1
    C.rho_power = 2.0
    C.omega_c = -15000
    
    C.SetHyperConstants(2)
    
    C.Csp = 380.
    
    C.phi_gamma = 3.0
    C.phi_z=3.5
    
    wr = Wrapper(C)
#     wr.solve(f0 = C.f0)
    C.SetHyperConstants(2)
    C.set_xs(np.array([0., 0., -28., 30., 30., 30., -15., -15.]))
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def waleckaMatsuiCut():
    C= eos.KVOR_mod2()
    C.Csp = 1.
    C.SetHyperConstants(2)
    C.Cs = 266.9
    C.Co = 195.7
    C.Cr = 54.71
    C.b = 0
    C.c = 0
    C.d = 0
    C.n0 = ((1.42*197.33/135)**3) / ((3*pi**2)/2.)
    C.f0 = 0.44
    
    C.alpha = 0
    C.omega_f = 0.5
    C.omega_kind = 1
    C.omega_a = 200.
    C.gamma = 0
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def myModFinal():
    C = myModR()
    C.rho_kind = 6
    C.rho_ld = 0
    C.rho_sat_a = 10
    C.rho_sat_f1 = 1
    C.rho_sat_f2 = 1
    C.gamma = 0.
    
    C.beta = 3.80530423
    C.rho_sat_val = 23.53502497
    C.rho_f = 0.522
    rho_f = 0.522
    C.rho_d = -4
    C.rho_width_f = C.f0
    C.rho_a = 0.479673
    C.rho_width_power = 1.
    C.rho_a0 = -2
    C.rho_a1 = 3
    C.rho_a2 = 0.8
    C.rho_a3 = -0
    C.rho_a4 = 0
    
    C.rho_e = 6
    
    C.beta1 = 0*1.8
    C.c1 = 80.
    
    C.omega_c = 0

    C.omega_kind = 2
    C.omega_a = 0.10669423
    C.omega_b = 7.09984178
    C.omega_f = 0.95
    wr = Wrapper(C)
    wr.solve(f0=C.f0, K0=240., J0=30.)
    return C

def myModNoPhi():
    C1 = myMod()
    wr1 = Wrapper(C1)
    C1.rho_ld = 0
    C1.rho_ld = 0
    C1.rho_sat_a = 0 
    C1.rho_sat_f1 = 1
    C1.rho_sat_f2 = 1
    C1.gamma = 0.
    
    
    C1.beta = 0
    C1.rho_sat_val = 0
    C1.rho_f = 0.522
    rho_f = 0.522
    C1.rho_d = 0
    C1.rho_width_f = C1.f0
    C1.rho_a = 1
    C1.rho_width_power = 0
    C1.rho_a0 = 0.5
    C1.rho_a1 = 2
    C1.rho_a2 = 4
    C1.rho_a3 = 0
    C1.rho_a4 = 1
    
    C1.rho_e = 0
    
    C1.rho = 0
    C1.drho = 0
    C1.d2rho = 0
    C1.d2om = 0
    C1.dom = 0
    C1.gamma2 = 0
    
    C1.rho_tan_a = 0
    C1.rho_tan_b = 0
    C1.rho_tan_c = 0
    
    C1.beta1 = 0*1.8
    C1.beta2 = 0
    C1.c1 = 80.
    

    C1.rho_tan_a = 3.

    C1.rho_a3 = 2
    C1.rho_a = 0.57682427
    C1.rho_sat_f1 = 0.43248055
    C1.rho_tan_a = 3.
    C1.rho_a0 = -1.06792546
    C1.rho_a1 = 9.72253398
    C1.rho_tan_b = 9.84217103
    C1.rho_a2 = 4
    C1.rho_a4 = 1
    C1.rho_tan_c = 0
    frange = np.linspace(0, 1, 100)
#     plt.plot(frange, map(C1.eta_r, frange))
#     plt.show()
    wr1.solve(f0 = C1.f0, K0=240., J0=30.)
    return C1           
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
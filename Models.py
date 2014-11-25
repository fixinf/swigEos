from Wrapper import Wrapper
import eosWrap as eos
import numpy as np

def myMod():
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
    C.b = 0
    C.c = 0
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def KVOR_cut():
    C = eos.KVOR_mod2()
    C.Csp = 1.
    C.Cs = 179.56233875171566
    C.Co =87.5996397368707
    C.Cr = 100.63642242798424
    C.b = 0.00773460805148428
    C.c = 0.00034461786646922604
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
    
    C.omega_a = 200.
    C.omega_f = 0.387
    
    C.phi_gamma = 3.
    C.phi_z = 3.5
    wr = Wrapper(C)
    wr.solve(f0=C.f0)
    C.SetHyperConstants(2)
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
    C.omega_f = 0.202
    
    C.phi_gamma = 3.
    C.phi_z = 3.5
    wr = Wrapper(C)
    wr.solve(f0=C.f0)
    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C


def KVOR_cut_02():
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
    C.omega_f = 0.202
    
    C.phi_gamma = 3.
    C.phi_z = 3.5
    wr = Wrapper(C)
    wr.solve(f0=C.f0)
    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def KVOR_cut_05():
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
    C.omega_f = 0.495
    
    C.phi_gamma = 3.
    C.phi_z = 3.5
    wr = Wrapper(C)
    wr.solve(f0=C.f0)
    C.SetHyperConstants(2)
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
    
import matplotlib.pyplot as plt
from Wrapper2 import Model
from Wrapper import Wrapper
import eosWrap as eos
import numpy as np
from scipy.misc.common import derivative
from scipy import optimize
from math import pi


def ret_fun(f):
    return 



def _KVOR():
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

def KVOR():
    return Model(_KVOR)

def _KVORcut04():
    C = eos.KVOR_cut()
    C.SetHyperConstants(2)
    C.omega_kind = 2
    C.Cs = 179.5623289954
    C.Co = 87.5996301756
    C.Cr = 100.6364192718
    C.b = 0.0077346088
    C.c = 0.0003446263
    C.Csp = 1.
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))

    C.omega_a = -0.5
    C.omega_f = 0.454
    C.omega_b = 55.76
    return C

def _KVORcut03():
    C = eos.KVOR_cut()
    C.SetHyperConstants(2)
    C.omega_kind = 2
    C.Cs = 179.5639888271
    C.Co = 87.5996147120
    C.Cr = 100.6364173848
    C.b = 0.0077354511
    C.c = 0.0003416056
    C.z = 0.65
    C.Csp = 1.
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))

    C.omega_a = -0.5
    C.omega_f = 0.365
    C.omega_b = 46.78
    # Wrapper(C).solve()
    return C

def KVORcut03():
    return Model(_KVORcut03)

def _KVORcut02():
    C = eos.KVOR_cut()
    C.SetHyperConstants(2)
    C.omega_kind = 2
    C.Cs = 184.0636560083
    C.Co = 87.5944165214
    C.Cr = 100.6364896138
    C.b = 0.0099027985
    C.c = -0.0073183493
    C.Csp = 1.
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))

    C.omega_a = -0.2
    C.omega_f = 0.249
    C.omega_b = 74.55
    # Wrapper(C).solve()
    return C

def KVORcut02():
    return Model(_KVORcut02)

def KVORcut04():
    return Model(_KVORcut04)

def _waleckaMatsui():
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

def _myMod():
    C = eos.KVOR_mod2()
    C.Cs = 234.1472066994
    C.Co = 134.8845385898
    C.Cr = 81.8421168107
    C.b = 0.0046749526
    C.c = -0.0029742081
    C.f0 = 0.27
    C.d = -0.5
    C.alpha = 0.4
    C.z = 0.65

    C.omega_kind = 2
    C.omega_a = 0.11
    C.omega_b = 7.1
    C.omega_f = 0.9

    C.rho_kind = 6

    C.rho_width_f = C.f0
    C.rho_width_power = 1.

    C.beta = 3.11
    C.rho_sat_val = 28.4
    C.rho_f = 0.522
    C.rho_a = 0.448
    C.rho_a0 = -0.614
    C.rho_a1 = 3.
    C.rho_a2 = 0.8
    C.rho_d = -4.
    C.rho_e = 6.

    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    # wr = Wrapper(C)
    # wr.solve(f0=C.f0, E0=-16., K0=240., J0=30.)
    return C

def _HyperTest():
    C = eos.KVOR()
    C.Csp = 1.
    C.SetHyperConstants(2)
    C.Cs = 164.462
    C.Co = 54.6041
    C.Cr = 121.690
    C.z = 0
    C.b = 0.0202832
    C.c = 0.0471633
    C.n0 = 0.17 * (197.33/135.)**3
    C.f0 = 0.15
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def __myMod_L(beta, gamma):
    C = _myMod()
    C.rho_kind = 1
    C.rho_a = 0.
    C.rho_power = 2.
    C.gamma = gamma
    C.beta = beta
    Wrapper(C).solve(f0=C.f0, K0=240.,J0=30.)
    return C

def myMod_L(beta, gamma):
    def _myMod_L():
        return __myMod_L(beta, gamma)
    _myMod_L.__name__ = '_myMod_Lb=%1.2f g=%1.2f'%(beta, gamma)
    return Model(_myMod_L)

def myMod():
    return Model(_myMod)

def waleckaMatsui():
    return Model(_waleckaMatsui)

def HyperTest():
    return Model(_HyperTest)

def __Cubero_cut(c):
    C = eos.KVOR_cut_sigma()
    # K=250., f0=0.2
    # C.Cs = 196.3428132661
    # C.Co = 90.7681870142
    # C.Cr = 88.7261140316
    # C.b = 0.0089455122
    # C.c = 0.0077076578

    # K =200., f0=0.2
    # C.Cs = 213.2607824762
    # C.Co = 90.7681524182
    # C.Cr = 88.7261190772
    # C.b = 0.0150061019
    # C.c = -0.0124942770

    #K=200., f0=0.19
    C.Cs = 205.0913029877
    C.Co = 84.4290026333
    C.Cr = 89.6408647791
    C.b = 0.0166215109
    C.c = -0.0114022252

    # K=200., f0=0.3
    # C.Cs = 273.0836626079
    # C.Co = 153.6303271470
    # C.Cr = 78.4149744148
    # C.b = 0.0057865796
    # C.c = -0.0072525183

    #K=200., f0=0.13
    C.Cs = 143.7534593933
    C.Co = 39.8340610093
    C.Cr = 95.5473126488
    C.b = 0.0102107661
    C.c = 0.2689864144

    C.z = 0
    C.b_sigma = 70.
    C.a_sigma = +1.
    C.c_sigma = c
    return C

def Cubero_cut(c, K0=None, f0=None, J0=None, suffix=None):
    def _Cubero_cut():
        return __Cubero_cut(c)
    _Cubero_cut.__name__ = '%s%1.2f'%(_Cubero_cut.__name__, c)
    return Model(_Cubero_cut, K0=K0, f0=f0, J0=J0, suffix=suffix)
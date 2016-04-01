import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from Wrapper2 import Model
from Wrapper import Wrapper
import eosWrap as eos
import numpy as np
from scipy.misc.common import derivative
from scipy import optimize
from math import pi
import Models

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

def _KVORd():
    C = eos.KVOR_d()
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

def KVOR_d():
    return Model(_KVORd, basefolder_suffix='KVOR_Delta')

def _KVORcut04():
    C = eos.KVORcut_d()
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
    C = eos.KVORcut_d()
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
    return Model(_KVORcut03, basefolder_suffix='KVORcut03_Delta')

def _KVORcut02():
    C = eos.KVORcut_d()
    # print([C.X_s[i] for i in range(12)])

    # C.SetHyperConstants(2)
    # print([C.X_s[i] for i in range(12)])

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

def _MS():
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

def _MKVOR():
    C = eos.MKVOR()
    C.Cs = 234.1472066994
    C.Co = 134.8845385898
    C.Cr = 81.8421168107
    C.b = 0.0046749526
    C.c = -0.0029742081
    C.f0 = 0.27
    C.d = -0.5
    C.alpha = 0.4
    C.z = 0.65

    C.a_om = 0.11
    C.b_om = 7.1
    C.f_om = 0.9

    C.beta = 3.11
    C.gamma = 28.4
    C.f_r = 0.522
    C.a_r0 = 0.448
    C.a_r1 = -0.614
    C.a_r2 = 3.
    C.a_r3 = 0.8
    C.d_r = -4.
    C.e_r = 6.

    C.SetHyperConstants(2)
    C.set_hs_alpha(np.array([0. for i in range(8)]))
    C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def MKVOR():
    return Model(_MKVOR)

def _MKVOR2():
    C = eos.MKVOR2()
    C.Cs = 234.1472066994
    C.Co = 134.8845385898
    C.Cr = 81.8421168107
    C.b = 0.0046749526
    C.c = -0.0029742081
    C.f0 = 0.27
    C.d = -0.5
    C.alpha = 0.4
    C.z = 0.65

    C.a_om = 0.11
    C.b_om = 7.1
    C.f_om = 0.9

    C.beta = 3.11
    C.gamma = 28.4
    C.f_r = 0.522
    C.a_r0 = 0.448
    C.a_r1 = -0.614
    C.a_r2 = 3.
    C.a_r3 = 0.8
    C.d_r = -4.
    C.e_r = 6.

    # C.SetHyperConstants(2)
    # C.set_hs_alpha(np.array([0. for i in range(8)]))
    # C.set_hs_z(np.array([0. for i in range(8)]))
    return C

def MKVOR2():
    M = Model(_MKVOR2)

    U = -50.
    xs_d = M.getDeltaXs(U)
    print(xs_d)
    M.setDeltaConst(np.array([xs_d for i in range(4)]),
             np.array([1. for i in range(4)]),
             np.array([1., 1., 1., 1.]),
             's = %.2f U = %.0f' % (xs_d, U))
    return M

def _MKVOR2_fom(f):
    C = _MKVOR2()
    C.fcut_om = f
    C.bcut_om = 100.

    return C

def MKVOR2_fom(f):
    __MKVOR2_fom = lambda: _MKVOR2_fom(f)
    __MKVOR2_fom.__name__ = _MKVOR2_fom.__name__ + '%.2f'%f
    M = Model(__MKVOR2_fom)
    U = -50.
    xs_d = M.getDeltaXs(U)
    print(xs_d)
    M.setDeltaConst(np.array([xs_d for i in range(4)]),
             np.array([1. for i in range(4)]),
             np.array([1., 1., 1., 1.]),
             's = %.2f U = %.0f' % (xs_d, U))
    return M

def _myModExp():
    return Models.myMod()

def MKVORexp():
    return Model(_myModExp)

def _MKVOR_d():
    C = eos.MKVOR_d()
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
    C.Cs = 234.1472110441
    C.Co = 134.8845431378
    C.Cr = 81.7485399957
    C.b = 0.0046749523
    C.c = -0.0029742083
    Wrapper(C).solve(f0=C.f0, K0=240.,J0=30.)
    return C

def _myModExp():
    C = _MKVOR_d()
    C.rho_kind = 1
    C.beta = 1.74991669
    C.gamma = 5.67486884
    C.rho_f = 0.42684265
    C.rho_a = 65.50171728
    return C

def myModExp():
    return Model(_myModExp)

def _myModExpOmega(omega_f):
    C = _myModExp()
    C.omega_kind = 3
    C.omega_a2 = 100.
    C.omega_f2 = omega_f
    return C

def myModExpOmega(omega_f):
    func = lambda: _myModExpOmega(omega_f)
    func.__name__ = _myModExpOmega.__name__ + '%.2f'%omega_f
    return Model(func)

def _MKValpha03(omega_f):
    C = _myModExpOmega(omega_f)
    C.alpha = 0.3
    C.Cs = 236.3405983248
    C.Co = 134.8845427865
    C.Cr = 81.7485395715
    C.b = 0.0047761455
    C.c = -0.0024441546
    Wrapper(C).solve(f0=C.f0, K0=240.,J0=30.)
    return C

def MKValpha03(omega_f):
    func = lambda: _MKValpha03(omega_f)
    func.__name__ = _MKValpha03.__name__ + '%.2f'%omega_f
    return Model(func)

def _MKValpha02(omega_f):
    C = _myModExpOmega(omega_f)
    C.alpha = 0.2
    C.Cs = 238.5606691067
    C.Co = 134.8845334593
    C.Cr = 81.7485399722
    C.b = 0.0048769734
    C.c = -0.0019268340
    Wrapper(C).solve(f0=C.f0, K0=240.,J0=30.)
    return C

def MKValpha02(omega_f):
    func = lambda: _MKValpha02(omega_f)
    func.__name__ = _MKValpha02.__name__ + '%.2f'%omega_f
    return Model(func)

def _MKValpha00(omega_f):
    C = _myModExpOmega(omega_f)
    C.alpha = 0.0
    C.Cs = 243.0809159724
    C.Co = 134.8845428070
    C.Cr = 81.7485403400
    C.b = 0.0050770283
    C.c = -0.0009275959
    # Wrapper(C).solve(f0=C.f0, K0=240.,J0=30.)
    return C

def MKValpha00(omega_f):
    func = lambda: _MKValpha00(omega_f)
    func.__name__ = _MKValpha00.__name__ + '%.2f'%omega_f
    return Model(func)

def myMod_L(beta, gamma):
    def _myMod_L():
        return __myMod_L(beta, gamma)
    _myMod_L.__name__ = '_myMod_Lb=%1.2f g=%1.2f'%(beta, gamma)
    return Model(_myMod_L)

def _OmegaWide(omega_f):
    C = _myModExpOmega(omega_f)
    C.omega_a2 = 30.
    return C

def OmegaWide(omega_f):
    func = lambda: _OmegaWide(omega_f)
    func.__name__ = _OmegaWide.__name__ + '%.2f'%omega_f
    return Model(func)

def myMod():
    return Model(_myMod)

def MKVOR_d():
    return Model(_MKVOR_d)

def waleckaMatsui():
    return Model(_waleckaMatsui)

def HyperTest():
    return Model(_HyperTest)

def __Cubero_cut(c, params=None):
    C = eos.KVORcut_sigma_d()
    # K=250., f0=0.2
    C.f0 = 0.2
    if params is not None:
        C.Cs = params['Cs']
        C.Co = params['Co']
        C.Cr = params['Cr']
        C.b = params['b']
        C.c = params['c']
    else:
        C.Cs = 196.3428132661
        C.Co = 90.7681870142
        C.Cr = 88.7261140316
        C.b = 0.0089455122
        C.c = 0.0077076578

    C.z = 0.
    C.b_sigma = 70.
    C.a_sigma = +1.
    C.c_sigma = c
    C.SetHyperConstants(2)
    return C

def __Cubero_cut2(c, params=None):
    C = eos.KVOR_cut_sigma2()
    # K=250., f0=0.2
    if params is not None:
        C.Cs = params['Cs']
        C.Co = params['Co']
        C.Cr = params['Cr']
        C.b = params['b']
        C.c = params['c']
    else:
        C.Cs = 196.3428132661
        C.Co = 90.7681870142
        C.Cr = 88.7261140316
        C.b = 0.0089455122
        C.c = 0.0077076578

    C.z = 0
    C.b_sigma = 11.1942
    C.a_sigma = 0.000269389 * (938 / 135)**4
    C.c_sigma = c
    C.SetHyperConstants(2)
    return C

def __Cubero_cut3(c, params=None):
    C = eos.KVOR_cut_sigma3()
    # K=250., f0=0.2
    if params is not None:
        C.Cs = params['Cs']
        C.Co = params['Co']
        C.Cr = params['Cr']
        C.b = params['b']
        C.c = params['c']
    else:
        C.Cs = 196.3428132661
        C.Co = 90.7681870142
        C.Cr = 88.7261140316
        C.b = 0.0089455122
        C.c = 0.0077076578

    C.z = 0
    C.b_sigma = 120.
    C.a_sigma = (139/135)**4
    C.c_sigma = c
    C.SetHyperConstants(2)
    return C

def __Cubero_cut_om(c, params=None):
    C = eos.KVOR_cut()
    C.a_sigma = 0
    # K=250., f0=0.2
    if params is not None:
        C.Cs = params['Cs']
        C.Co = params['Co']
        C.Cr = params['Cr']
        C.b = params['b']
        C.c = params['c']
    else:
        C.Cs = 196.3428132661
        C.Co = 90.7681870142
        C.Cr = 88.7261140316
        C.b = 0.0089455122
        C.c = 0.0077076578

    C.z = 0
    C.omega_b = 70.
    C.omega_a = -0.5
    C.alpha = 0
    C.omega_f = C.f0 + (1 - C.f0) * c
    # print(C.f0 + (1-C.f0)*C.c_omega, C.a_omega, C.b_omega)
    C.SetHyperConstants(2)
    return C


def Cubero_cut_om(c, K0=None, f0=None, J0=None, suffix=None, params=None):
    def _Cubero_cut_om():
        return __Cubero_cut_om(c, params)
    _Cubero_cut_om.__name__ = _Cubero_cut_om.__name__ + '%.2f'%c + suffix
    return Model(_Cubero_cut_om, K0=K0, f0=f0, J0=J0, suffix=suffix, basefolder_suffix='NLWCutsOm/K%2.fm%1.1f'%(K0, 1-f0))

def Cubero_cut(c, K0=None, f0=None, J0=None, suffix=None, params=None):
    def _Cubero_cut():
        return __Cubero_cut(c, params)
    _Cubero_cut.__name__ = _Cubero_cut.__name__ + '%.2f'%c + suffix
    return Model(_Cubero_cut, K0=K0, f0=f0, J0=J0, suffix=suffix, basefolder_suffix='NLWCuts_Delta')

def Cubero_cut2(c, K0=None, f0=None, J0=None, suffix=None, params=None):
    def _Cubero_cut():
        return __Cubero_cut2(c, params)
    _Cubero_cut.__name__ = _Cubero_cut.__name__ + '%.2f'%c + suffix
    return Model(_Cubero_cut, K0=K0, f0=f0, J0=J0, suffix=suffix, basefolder_suffix='NLWCuts2/K%2.fm%1.1f'%(K0, 1-f0))

def Cubero_cut3(c, K0=None, f0=None, J0=None, suffix=None, params=None):
    def _Cubero_cut():
        return __Cubero_cut3(c, params)
    _Cubero_cut.__name__ = _Cubero_cut.__name__ + '%.2f'%c + suffix
    return Model(_Cubero_cut, K0=K0, f0=f0, J0=J0, suffix=suffix, basefolder_suffix='NLWCuts3/K%2.fm%1.1f'%(K0, 1-f0))


def __Cubero_cut_rho(c, params=None):
    C = eos.KVOR_cut_rho()
    C.a_sigma = 0
    # K=250., f0=0.2
    if params is not None:
        C.Cs = params['Cs']
        C.Co = params['Co']
        C.Cr = params['Cr']
        C.b = params['b']
        C.c = params['c']
    else:
        C.Cs = 196.3428132661
        C.Co = 90.7681870142
        C.Cr = 88.7261140316
        C.b = 0.0089455122
        C.c = 0.0077076578

    C.z = 0
    C.b_rho = 70.
    C.rho_a = -0.5
    C.alpha = 0
    C.rho_f = C.f0 + (1 - C.f0) * c
    # print(C.f0 + (1-C.f0)*C.c_omega, C.a_omega, C.b_omega)
    C.SetHyperConstants(2)
    return C


def myMod2():
    def _myMod2():
        return _myMod()
    return Model(_myMod2)


def Cubero_cut_rho(c, K0=None, f0=None, J0=None, suffix=None, params=None):
    def _Cubero_cut_rho():
        return __Cubero_cut_rho(c, params)
    _Cubero_cut_rho.__name__ = _Cubero_cut_rho.__name__ + '%.2f'%c + suffix
    return Model(_Cubero_cut_rho, K0=K0, f0=f0, J0=J0, suffix=suffix, basefolder_suffix='NLWCutsRho/K%2.fm%1.1f'%(K0, 1-f0))


def _Wal_d():
    C = eos.Walecka_d()
    C.Cs = 196.3428132661
    C.Co = 90.7681870142
    C.Cr = 88.7261140316
    C.b = 0.0089455122
    C.c = 0.0077076578
    C.n0 = 0.16*(197.33/135)**3
    C.f0 = 0.2
    return C

def Wal_d(J0=None, f0=0.2):
    return Model(_Wal_d, K0=250., J0=J0, f0=f0, basefolder_suffix='Walecka_Delta')

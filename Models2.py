import matplotlib.pyplot as plt
from Wrapper2 import Model
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

def waleckaMatsui():
    return Model(_waleckaMatsui)



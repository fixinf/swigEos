import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join

models = [Models.KVOR, Models.myMod]
for model in models:
    C = model()
    
    wr = Wrapper(C)
    n, res = wr.dumpJ()
    
    f0prime = wr.f0prime(n)
    pf = np.array(map(eos.p_f, n/2))
    Es, fs = wr.Esymm(n, ret_f=1)
    Ef = np.sqrt(pf**2 + C.M[0]**2 * np.array(map(lambda z: C.phi_n(0, z)**2, fs)))
    a4 = wr.m_pi * pf**2 / (6 * Ef) * (1 + f0prime)
    
    plt.plot(n/wr.n0, a4)
#     plt.plot(n/wr.n0, res)
plt.show()


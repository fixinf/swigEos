import eosWrap as eos
import numpy as np
from matplotlib import pyplot as plt

import Models2

C = Models2.KVOR().C

init = np.array([0.])
finit = np.array([0.])

n = C.n0

res = eos.stepE_rho(n, init, finit, len(finit), 30, 0., C)
res2 = eos.stepE(n, init, finit, len(finit), 30, C)

print(res)
print(res2)
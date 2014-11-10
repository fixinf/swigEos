import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt


C = eos.KVOR()
C.SetHyperConstants(2)
wr = Wrapper(C)
C.Cs = 200.
wr.solve()
wr.reset(nmin=0., nmax = 4., npoints = 2000, iter=100)

wr.setDriver()

N, M, R = wr.stars(nmin=0.5, nmax=4., npoints=100)
with open('KVORMasses.dat', 'w') as f:
    for i, _n in enumerate(N):
        f.write('%f %f %f \n'%(_n/wr.n0, M[i], R[i]))
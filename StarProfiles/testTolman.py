import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join


model = Models.KVOR
C = model()
wr = Wrapper(C)
N, M, R, E, P = wr.dumpStarProfile(2.5*wr.n0)
# plt.plot(N, E)
tolmanE = E[0] - 5*(P[0] - P) + 8 * (P[0] - P)**2 / (E[0] + P[0])
plt.plot(P, E)
plt.show() 
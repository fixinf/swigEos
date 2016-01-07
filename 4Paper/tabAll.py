import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models2
from tabulate import tabulate
from os.path import join
import inspect

# models = [Models.myModChi]#, Models.KVOR_cut_03Chi, Models.KVORphi]
# models = [Models.KVORphi]
# models = [Models.KVOR_cut_03Chi]

models = [Models2.KVOR(), Models2.KVORcut03(), Models2.myMod()]
DUs = []

for wr in models:
    wr.dumpEos()
    # wr.hyper_phi.dumpEos()
    # wr.hyper_phi_sigma.dumpEos()

exit()

import eosWrap as eos
import numpy as np
from matplotlib import pyplot as plt

import Models2

wr = Models2.MKValpha00(0.66)
wr.dumpAll(hyper=0)
exit()
m = wr.delta_phi
m.dumpDeltaSym()
m.dumpEos()
m.dumpMassesCrust()

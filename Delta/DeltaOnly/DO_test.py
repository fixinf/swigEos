import eosWrap as eos
import Models2
import numpy as np
from matplotlib import pyplot as plt

wr = Models2.KVOR_d()

m = wr.delta_only

m.dumpEos()

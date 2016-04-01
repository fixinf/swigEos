import eosWrap as eos
import numpy as np
from matplotlib import pyplot as plt
import Models2
from scipy.optimize.minpack import curve_fit

for f in np.arange(0.72, 0.8, 0.04):
    wr = Models2.MKVOR2_fom(f)
    wr.dumpProps()
    wr.dumpScalings()
    wr.dumpEos()
    wr.nucl.dumpEos()
    # plt.plot(wr.nrange, wr.nucl._P)
    # plt.show()
    # plt.plot(wr.nrange, wr.nucl._E)
    # plt.show()
    # wr.dumpAll(hyper=0)
    wr.delta_sym.dumpEos()


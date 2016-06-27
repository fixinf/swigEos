import Models2
import numpy as np

from matplotlib import pyplot as plt
# wr = Models2.MKVOR_d()
wr = Models2.MKVOR2final()
# for m in [wr.rcond_nucl, wr.rcond_hyper_phi, wr.rcond_hyper_phi_sigma,
#           wr.rcond_delta_phi, wr.rcond_delta_phi_sigma]:

for m in [wr.rcond_delta_phi_sigma]:
      m.dumpEos()
      m.dumpMassesCrust() 
exit()

for wr in [wr1, wr2, wr3]:
    # wr.dumpAll(hyper=0)
    wr.delta_sym.dumpEos()
    # wr.delta_only.dumpEos()
    # wr.delta_only.dumpMassesCrust()
    # wr.dumpBaryonParams()
    # exit()
    # U = -50.
    # xs = wr.getDeltaXs(U)

    # wr.setDeltaConst(np.array([xs for i in range(4)]), 
    #     np.array([1. for i in range(4)]),
    #     np.array([1. for i in range(4)]),
    #     's = %.2f U = %.0f' % (xs, U))
    # exit()
    # wr.nucl.dumpEos()
    # wr.nucl.loadEos()
    # wr.nucl.dumpMassesCrust()
    # plt.plot(wr.nucl._E, wr.nucl._P)
    # plt.show()
    # exit()

    # wr.nucl.dumpEos()
    # wr.nucl.dumpMassesCrust()
    # exit()
    # wr.dumpAll(hyper=0)
    # for m in [wr.delta_phi, wr.delta_phi_sigma]:
    #     print(m.foldername)
    #     # exit()
    #     m.dumpEos()

    #     # m.dumpMassesCrust()

    # # exit()

    # wr.dumpAll(hyper=0)


    # for m in [wr.nucl, wr.hyper_phi, wr.hyper_phi_sigma]:
    # 	m.dumpEos()
    # 	# m.dumpMassesCrust()
    # # exit()
    # for m in [wr.rcond_nucl, wr.rcond_hyper_phi, wr.rcond_hyper_phi_sigma]:
    # 	m.dumpEos()
    # 	# m.dumpMassesCrust()
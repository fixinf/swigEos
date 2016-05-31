import Models2
import numpy as np
wr = Models2.KVORcut03()


xs = wr.getDeltaXs(-100.)
xr = 1.
xo = 1.
wr.setDeltaConst(np.array([xs, xs, xs, xs]),
                 np.array([1., 1., 1., 1.]),
                 np.array([xr, xr, xr, xr]),
                 "_xs=%.2f_xo=%.2f_xr=%.2f" % (xs, xo, xr))
wr.delta_sym.dumpEos()
# wr.delta_phi.dumpEos()
# wr.delta_phi.dumpMassesCrust()
# wr.delta_phi_sigma.dumpEos()
# wr.delta_phi_sigma.dumpMassesCrust()

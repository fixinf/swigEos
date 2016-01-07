__author__ = 'const'
import Models2
import numpy as np
# m = Models2.KVORcut03()
m = Models2.MKVOR_d()
# m = Models2.Wal_d()
for xs, xo in [[1., 1.]]:
    m.setDeltaConst(np.array([xs for i in range(4)]),
                     np.array([xo for i in range(4)]),
                     np.array([1., 1., 1., 1.]),
                     's=%.2f o=%.2f'%(xs, xo))
    # m.delta_sym.dumpEos()
    # m.delta.dumpEos()
    # m.delta_phi.dumpEos()
    m.delta_phi_sigma.dumpEos()
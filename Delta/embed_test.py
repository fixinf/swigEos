__author__ = 'const'
import Models2
import numpy as np
m = Models2.MKVOR_d()

xs = 1.31
xo = 1.0

m.setDeltaConst(np.array([xs for i in range(4)]),
                 np.array([xo for i in range(4)]),
                 np.array([1., 1., 1., 1.]),
                 's=%.2f o=%.2f'%(xs, xo))

s = m.delta_sym
print(s.foldername)
print(s.n0, s.C.n0)
# exit()
s.dumpEos()
s.dumpNc()
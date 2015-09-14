from scipy import optimize

__author__ = 'const'
import Models2
import numpy as np
from matplotlib import pyplot as plt

wr = Models2.KVOR_d()

m = wr.delta

xs = 1.25
xo = 1.

for xs, xo in [[1.25, 1.], [1., 1.], [1.15, 0.9]]:
    wr.setDeltaConst(np.array([xs for i in range(4)]),
                     np.array([xo for i in range(4)]),
                     np.array([1., 1., 1., 1.]),
                     's=%.2f o=%.2f'%(xs, xo))

    m.dumpDeltaSym('Ebind')
exit()

# def n_crit(xs, xo, n_init=0.):
#     wr.setDeltaConst(np.array([xs for i in range(4)]),
#                      np.array([xo for i in range(4)]),
#                      np.array([1., 1., 1., 1.]),
#                      '1')
#
#     f0 = 0.195
#
#     # print([m.C.X_s[i] for i in range(12)])
#     _n = 0.
#     for n in np.linspace(m.nrange[1], m.nrange[-1], 100):
#         if m.getSymm(n)[0] > 1e-6:
#             print(n)
#             return n
#
# # print(n_crit(1.4, .9))
#
# # x_range = np.linspace(0, 2, 100)
# # nclist = []
# # for x in x_range:
# #
# # plt.plot(x_range, n_c)
# # plt.show()
# # exit()

res = []
lastn = 0.
lastf = 0.



for i, n in enumerate(m.nrange):
    _res = m.getSymm(n, lastx=lastn, lastf=lastf)
    print('res = ', _res)
    if (i > 1):
        lastn = _res[0][0]
    print('lastn =', lastn)
    lastf = _res[-1]
    res.append(_res)

res = np.array(res)

n_d = res[:, 0]
f = res[:, 2]

print(n_d)
print(f)

plt.plot(m.nrange/m.n0, n_d/m.nrange)
plt.plot(m.nrange/m.n0, f)
plt.ylim([0,3])
plt.show()


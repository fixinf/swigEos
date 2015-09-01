from os.path import join
from scipy import optimize
from scipy.interpolate.interpolate import interp1d
from matplotlib import pyplot as plt
from scipy.optimize.zeros import brentq, bisect

__author__ = 'const'
import Models2
import numpy as np

model = Models2.KVORcut03()
for m in [model.hyper, model.hyper_phi, model.hyper_phi_sigma]:
    # m = Models2.KVOR().hyper_phi
    inter_kind = 'cubic'
    data = np.loadtxt(join(m.foldername, m.filenames['mass_crust']), skiprows=1)
    n = data[:, 0]
    iData = [interp1d(n, data[:, i], kind=inter_kind) for i in range(1, data.shape[1])]
    iStr = [interp1d(n, data[:, i], kind=inter_kind) for i in range(6, data.shape[1])]
    maxRes = optimize.minimize_scalar(lambda z: -iData[0](z), bounds=(0.1, 9), method='bounded', tol=1e-16)

    mMax = -maxRes.fun
    nMax = maxRes.x

    # plt.plot(n, iData[0](n))
    # fig = plt.gcf()
    # fig.gca().add_artist(plt.Circle((nMax, mMax), radius=.01))
    # plt.show()

    print(m.__class__)
    print(nMax, mMax)
    # strMax = (iStr[0](mMax) + iStr[1](mMax) + iStr[2](mMax) + iStr[3](mMax) + 2 * iStr[4](mMax) + 2*iStr[5](mMax))
    strMax = (iStr[0](nMax) + iStr[1](nMax) + iStr[2](nMax) + iStr[3](nMax) + 2 * iStr[4](nMax) + 2*iStr[5](nMax))
    print(strMax)
    print('DU thresholds:')
    nDu = []
    for s in iStr:
        try:
            _ndu = bisect(lambda z: s(z) - 1e-6, 1, 8, xtol=1e-12




            )
        except ValueError:
            _ndu = 0
        nDu.append(_ndu)
    print(nDu)
    print('DU masses:')
    print([iData[0](_n) if _n > 1e-6 else 0 for _n in nDu])
    print('\n')
# plt.plot(data[:, 1], np.array([istr(data[:, 1]) for istr in iStr]).transpose())
# plt.show()



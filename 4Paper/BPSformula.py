import numpy as np
from matplotlib import pyplot as plt
import Models2


wr = Models2.KVOR()


### 1 MeV/fm^3 = 1.7827e12 g/cm^3
mevfm3 = 1.7827e12
### MeV in g
mev = 1.7827e-27 #g
### Nucleon mass in HP: m0 = 1.66e-24 g
m0 = 1.66e-24

def f0(x):
    return 1./(1+ np.exp(x))

def getEN(n):
    """
    :param n: baryon number density in fm^-3
    :return: energy density for the HP interpolation formula in MeV/fm^3
    """
    p = [
        0.320,
        2.17,
        0.173,
        3.01,
        0.540,
        0.847,
        3.581
    ]

    q = [
        0.608,
        2.41,
        2.39,
        3.581,
        1.681,
        0.850,
        11.64
    ]

    # p = [
    #     0.423,
    #     2.42,
    #     0.031,
    #     0.78,
    #     0.238,
    #     0.912,
    #     3.674
    # ]
    #
    # q = [
    #     0.183,
    #     1.26,
    #     6.88,
    #     3.612,
    #     2.248,
    #     0.911,
    #     11.56
    # ]

    res = (1. + (p[0] * n**p[1] + p[2] * n**p[3])/(1 + p[4]*n)**2 * f0(-p[5] * (np.log10(n) + p[6])) +
           n / (8e-6 + 2.1 * n ** 0.585) * f0(p[5] * (np.log10(n) + p[6])))

    print(m0/mev)
    res *= m0*n / mev # mev/fm**3
    return res

Ec, Pc, Nc = np.loadtxt('ex_bps_input.txt', skiprows=1).transpose()
Nc /= 1e39
fps = np.loadtxt('testfps.dat', skiprows=0)

nrange = np.linspace(0., 0.16, 100)
print(getEN(0.1581849930E-08) * mev * 1e39)
print(getEN(0.3977277962) * mev * 1e39)
# exit()
print(fps)
plt.plot(fps[:, 1], fps[:, 2], marker='o')
plt.plot(nrange, getEN(nrange) * mev * 1e39)
plt.plot(Nc, Ec, marker='o')
plt.xlim([0., 0.17])
plt.ylim([0., 0.1e16])

plt.show()

plt.semilogx(fps[:, 1], fps[:, 2], marker='o')
plt.semilogx(nrange, getEN(nrange) * mev * 1e39)
plt.semilogx(Nc, Ec, marker='o')
plt.xlim([0., 0.17])
plt.ylim([0., 0.1e16])

plt.show()

plt.plot(fps[:, 1], fps[:, 2]/fps[:, 1], marker='o')
plt.plot(nrange, getEN(nrange) * mev * 1e39 / nrange)
plt.plot(Nc, Ec / Nc, marker='o')

plt.show()

plt.semilogx(fps[:, 1], fps[:, 2]/fps[:, 1], marker='o')
plt.semilogx(nrange, getEN(nrange) * mev * 1e39 / nrange)
plt.semilogx(Nc, Ec / Nc, marker='o')

plt.show()



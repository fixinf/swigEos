__author__ = 'const'
import matplotlib
matplotlib.use('QT4Agg')
from matplotlib import pyplot as plt
from cmath import log
from math import sqrt, pi
import numpy as np
import Models2
import Models
from Wrapper import Wrapper
import eosWrap as eos


def I1(m, x):
    t = sqrt(1 + x**2)
    return 2 * m**2 * (0.5*x*t + x/t - 1.5 * log(x + t))/pi**2


def I2(m, x):
    t = sqrt(1 + x**2)
    return 0.5 * m**4 * (x*t**3 - 0.5*x*t - 0.5 * log(x + t))/pi**2


def I3(m, x):
    t = sqrt(1 + x**2)
    return m**3 * (x*t - log(x + t))/pi**2


def eq(x, params):
    n0, E, m, meff, K, J = params
    Cs, b, c = x

    k = (3*pi**2 * n0/2)**(1/3)
    ek = sqrt(meff**2 + k**2)

    # print(Co, (m + E - ek)/n0)
    Co = (m + E - ek)/n0
    # print('Co: ',Co)

    alpha = [K - Co*6*k**3/pi**2 - 3*k**2/ek,
             0.5*(m - meff)**2,
             m - meff]

    beta = [2*(m-meff)*m*alpha[0],
            m*(m-meff)**3 / 3,
            m*(m-meff)**2]

    gamma = [3*(m-meff)**2 * alpha[0],
            0.25*(m - meff)**4,
            (m-meff)**3]

    delta = [-alpha[0]*I1(meff, k/meff) - 6 * k**3 * (meff/ek)**2/pi**2,
             n0*(m + E) - I2(meff, k/meff) - 0.5 * Co * n0**2,
             I3(meff, k/meff)]

    # print('111 ', Cs, alpha[0]/(delta[0] - gamma[0]*c - beta[0] * b))

    res = np.array([alpha, beta, gamma]).transpose()
    # print(alpha, beta, gamma)
    # print(res)
    # print(x)
    x[0]=1/x[0]
    # print(x)
    # print(delta)
    # print('eq1: ', res.dot(np.array(x)) - np.real(np.array(delta)))
    x1 = np.linalg.solve(res, np.real(np.array(delta)))
    # print(x1)
    # print('eq2: ', res.dot(np.array(x1)) - np.real(np.array(delta)))

    Cr = 8*(J - k**2/6/ek)/n0

    return list(x1) + [(m + E - ek)/n0] + [Cr]

#
#
# C = Models.Cubero2()
# wr = Wrapper(C)
# print(C.M[0]*135)
# mn = C.M[0]
# Cs = 12.684 * (1/197.33)**2 * 135**2
# Co = 7.148 * (135./197.33)**2
# b = 0.5610e-2
# c = -0.6986e-2
# E = -16.3
# n0 = 0.153/0.15 * C.n0
# K = 200
# f = 0.3
#
# # Cs = C.Cs/mn**2
# # Co = C.Co/mn**2
# # b = C.b
# # c = C.c
# # E = wr.ESbind(np.linspace(0, C.n0, 10))[-1]
# # print(E)
# # n0 = C.n0
# # K = wr.K()
# # print(K)
# # f = C.f0
#
# Csm1, b, c, Com1 = eq([Cs, b, c],
#    [Co, n0, E/135, C.M[0], C.M[0]*(1-f), K/135])
#
# # 197.33 Mev fm = 1
# # 1/fm = 197.33 MeV = 197.33/135 mpi
# # fm**2 = (135/197.33) mpi
#
# # C.Cs = (1/Csm1)*C.M[0]**2
# # C.b = b
# # C.c = c
# # C.Co = Com1*C.M[0]**2
# # print(C.Cs, C.Co)1
# wr.n0 = n0
#
# C.Cs = Csm1**-1 * mn**2
# C.Co = Com1 * mn**2
# C.b = b
# C.c = c
# C.f0 = f
# # C.f0 = 0.3
# E, flist = wr.Esymm(np.linspace(0, n0, 10), ret_f=1)
# print(flist)
# print(n0, np.linspace(0, n0, 10)[-1])
# print(wr.K(), wr.ESbind(np.linspace(0, n0, 5))[-1])
# nrange = np.linspace(0, 2*n0, 1000)
# plt.plot(nrange/n0, wr.ESbind(nrange))
# plt.show()
# # wr.solve(f0=0.3, E0=-16.3, K0=200., J0=32.5)










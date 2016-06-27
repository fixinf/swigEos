import Models2
import numpy as np
import joblib as jl
from matplotlib import pyplot as plt
from os.path import join
import multiprocess as mp

from pprint import pprint as print

from matplotlib import pyplot as plt

wr = Models2.MKVOR2final()

# wr.delta_sym.dumpNc(ulist)

# exit()~

def fun_sym(U):
    xs_d = wr.getDeltaXs(U)
    wr.setDeltaConst(np.array([xs_d for i in range(4)]),
                     np.array([1. for i in range(4)]),
                     np.array([1., 1., 1., 1.]),
                     'xs=%.2f U = %.2f' % (xs_d, U))
    print(wr.delta_phi.foldername)
    # exit()
    for m in [wr.delta_sym]:
        m.dumpEos()
        m.dumpMu()
        # m.dumpMassesCrust()
    return 0.


def fun(U):
    xs_d = wr.getDeltaXs(U)
    wr.setDeltaConst(np.array([xs_d for i in range(4)]),
                     np.array([1. for i in range(4)]),
                     np.array([1., 1., 1., 1.]),
                     'xs=%.2f U = %.2f' % (xs_d, U))

    print(wr.delta_phi.foldername)
    # exit()
    for m in [wr.delta_phi, wr.delta_phi_sigma]:
    # for m in [wr.delta_only]:
        m.dumpEos()
        m.dumpScalingsN()
        # m.dumpMassesCrust()
    return 0

def get(U, get_mass=1):
    xs_d = wr.getDeltaXs(U)
    wr.setDeltaConst(np.array([xs_d for i in range(4)]),
             np.array([1. for i in range(4)]),
             np.array([1., 1., 1., 1.]),
             'xs=%.2f U = %.2f'%(xs_d, U))
    print(wr.delta_phi.foldername)
    m = wr.delta_phi
    m.loadEos()
    if get_mass:
        mass = np.loadtxt(join(m.foldername, m.filenames['mass_crust']+'_linear'), skiprows=1)
        mmax = max(mass[:, 11])
    lept = m.lepton_concentrations(ret_du=1)
    conc = m.concentrations()
    # plt.plot(m.nrange/m.n0, lept[2])
    # plt.plot(m.nrange/m.n0, conc[:, 1])
    # plt.show()
    # plt.plot(m.nrange/m.n0, np.abs(lept[2] - conc[:, 1]))
    # plt.show()
    # print(np.argmin(np.abs(lept[2] - conc[:, 1])))
    n_du = m.nrange[np.argmin(np.abs(lept[2] - conc[:, 1]))]/m.n0
    _set = [0 for i in range(m.n_baryon - 2)]
    nc = [max(m.nrange/m.n0) for s in _set]
    for i, r in enumerate(m.concentrations()):
        # print(r)
        for j, s in enumerate(_set):
            if not s:
                if r[-len(_set)+j] > 1e-7:
                    nc[j] = m.nrange[i] / m.n0
                    _set[j] = 1
    if get_mass:
        return [mmax, n_du] + nc
    else:
        return [n_du] + nc


def get_sym(U):
    xs_d = wr.getDeltaXs(U)
    wr.setDeltaConst(np.array([xs_d for i in range(4)]),
                     np.array([1. for i in range(4)]),
                     np.array([1., 1., 1., 1.]),
                     'xs=%.2f U = %.2f' % (xs_d, U))
    print(wr.delta_sym.foldername)
    m = wr.delta_sym
    m.loadEos()
    # mass = np.loadtxt(join(m.foldername, m.filenames['mass_crust']+'_linear'), skiprows=1)
    # mmax = max(mass[:, 1])
    # lept = m.lepton_concentrations(ret_du=1)
    # plt.plot(m.nrange/m.n0, m.rho)
    conc = m.concentrations()
    # plt.plot(m.nrange/m.n0, lept[2])
    # plt.plot(m.nrange/m.n0, conc[:, 1])
    # plt.show()
    # plt.plot(m.nrange/m.n0, conc)
    # plt.show()
    nc = 8.
    _set = 0
    nm = 8.
    for i, r in enumerate(m.concentrations()):
        if not _set:
            if r[10] > 1e-7:
                nc = m.nrange[i] / m.n0
                _set = 1

    for i, r in enumerate(m.rho):
        # print(r)
        if not (r[0] < 1.):
            nm = m.nrange[i]/m.n0
            break
    return [nc, nm]

ulist = np.linspace(-50., -150., 20, endpoint=0)
# ulist = np.array([-50., -60., -70.,-90., -100., -110.,-130.,-150.,-157.])
# ulist = [-55., -60., -65]
# ulist = np.array([-50., -55., -60., -65., -70., -75., -80., -85., -90., -100., -110., -120., -130., -140., -150.,-157.])
# ulist = np.linspace(-50., -150., 50, endoint=0)
res = []
# ulist=[-150.]
# for MKVOR
# ulist = np.arange(-50., -131., -5)
for u in ulist:
# for u in [-148.]:
    # fun_sym(u)
    fun(u)
    # get_sym(u)
    # if u > -105:
    # res.append([u] + get_sym(u))
    # res.append([u] + get(u))
    # res.append([u] + get(u, get_mass=0))        
    # else:
    #     break
exit()
res = np.array(res)
print(res)
# np.savetxt(join(wr.foldername, 'masses_du_do.dat'), res)
print(join(wr.foldername, 'nc_dp_U.dat'))
np.savetxt(join(wr.foldername, 'nc_dp_U.dat'), res, fmt='%.6f')
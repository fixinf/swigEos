import Models2
import numpy as np
import joblib as jl
from matplotlib import pyplot as plt
from os.path import join
import multiprocess as mp
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from pprint import pprint as print

from matplotlib import pyplot as plt

# wr = Models2.KVORcut03()

# wr = Models2.MKVOR_d()
wr = Models2.MKVOR2_exp()

# wr = Models2.MKVOR2final()

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
        # m.dumpEos()
        # m.dumpScalingsN()

        m.dumpEos()
        if not m.needsMaxw():
            pass
            # m.dumpMassesCrust()
        else:
            print(m.foldername)
            # m.processMaxw()
            # plt.plot(m.nrange/m.n0, m._P)
            # plt.show()
            # m.dumpMassesCrust()
    return 0

def get(U, get_mass=1):
    xs_d = wr.getDeltaXs(U)
    wr.setDeltaConst(np.array([xs_d for i in range(4)]),
             np.array([1. for i in range(4)]),
             np.array([1., 1., 1., 1.]),
             'xs=%.2f U = %.2f'%(xs_d, U))
    print(wr.delta_phi.foldername)
    m = wr.delta_only

    # m = wr.delta_only
    m.loadEos()
    if get_mass:
        mass = np.loadtxt(join(m.foldername, m.filenames['mass_crust']+'_linear'), skiprows=1)
        mmax = max(mass[:, 1])
        iM = interp1d(mass[:, 0], -mass[:, 1], kind='cubic')
        nmmax = mass[np.argmax(mass[:, 1]), 0]
        nmmax = minimize(iM, nmmax).x
    else:
        mmax = 0.
        nmmax = 0.
    lept = m.lepton_concentrations(ret_du=1)
    conc = m.concentrations()
    # plt.plot(m.nrange/m.n0, lept[2])
    # plt.plot(m.nrange/m.n0, conc[:, 1])
    # plt.show()
    # plt.plot(m.nrange/m.n0, np.abs(lept[2] - conc[:, 1]))
    # plt.show()
    # print(np.argmin(np.abs(lept[2] - conc[:, 1])))

    # n_du = m.nrange[np.argmin(np.abs(lept[2] - conc[:, 1]))]/m.n0
    n_du, m_du = m.getDuCrit()
    n_du /= m.n0
    _set = [0 for i in range(m.n_baryon - 2)]
    nc = [max(m.nrange/m.n0) for s in _set]
    for i, r in enumerate(m.concentrations()):
        # print(r)
        for j, s in enumerate(_set):
            if not s:
                if r[-len(_set)+j] > 1e-7:
                    nc[j] = m.nrange[i] / m.n0
                    _set[j] = 1

    return [mmax, n_du] + nc + [m_du] + [nmmax]



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
    nc = max(wr.nrange)/wr.n0
    _set = 0
    nm = max(wr.nrange)/wr.n0
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
    fjump = 1.
    njump = 8.
    f = m.rho[:, 0]
    if any(np.diff(f) > 0.02):
        fjump = f[np.argmax(np.diff(f)) + 1]
        njump = m.nrange[np.argmax(np.diff(f)) + 1]/m.n0
    return [nc, nm, fjump, njump]

# ulist = np.linspace(-50., -150., 50, endpoint=0)
# ulist = np.arange(-130, -151, -5)
# ulist = np.arange(-146., -151., -2)
ulist = np.arange(-146., -151., -2)
# ulist = np.linspace(-60, -55, 10)
# ulist = np.arange(-100., -131., -5)
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
# print(join(wr.foldername, 'nc_dp_U.dat'))
# np.savetxt(join(wr.foldername, 'nc_do_U.dat'), res, fmt='%.6f')
np.savetxt(join(wr.foldername, 'nc_sym_U.dat'), res, fmt='%.6f')
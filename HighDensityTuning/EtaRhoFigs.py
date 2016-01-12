from os.path import join
from scipy.interpolate.interpolate import interp1d

import Models2
from matplotlib import pyplot as plt
import numpy as np

wr = Models2.myModExp()
m = wr.hyper
m.reset(timeout=6)
m.dumpEos()
n_V = np.array([
    sum([m.C.X_o[i]*n[i] for i in range(m.n_baryon)])
    for n in m.concentrations()
])

n_I = np.array([
    sum([m.C.X_r[i]*m.C.T[i]*n[i] for i in range(m.n_baryon)])
    for n in m.concentrations()
])

plt.plot(m.nrange[:len(n_I)]/m.n0, n_I)
plt.show()

iNv = interp1d(m.nrange[: len(n_V)], n_V)
iNi = interp1d(m.nrange[: len(n_I)], n_I)

flist = m.rho[:, 0]
Ueff = np.array([
    m.C.U(f) + m.C.Co*n_V[i]**2 / (2 * m.C.M[0]**2 * m.C.eta_o(f)) +
    m.C.Cr*n_I[i]**2 / (2 * m.C.M[0]**2 * m.C.eta_r(f))
    for i, f in enumerate(flist)
])

n = 2*m.n0
for n in [m.n0, 2*m.n0, 3*m.n0, 3.8*m.n0, 4*m.n0, 5*m.n0, 5.4*m.n0]:
    frange, eq = m.get_feq(n)
    np.savetxt(join(m.foldername, m.name+'_eq_f_%.1f.dat'%(n/m.n0)),
               np.array([frange, eq]).transpose())
    plt.plot(frange, eq)
plt.show()

Ueff_f = np.array([
    m.C.U(f) + m.C.Co*n_V[i]**2 / (2 * m.C.M[0]**2 * m.C.eta_o(f)) +
    m.C.Cr*n_I[i]**2 / (2 * m.C.M[0]**2 * m.C.eta_r(f))
    for i, f in enumerate(frange)
])

Ueff_list = np.array([
    [
        m.C.U(f) + m.C.Co*nv**2 / (2 * m.C.M[0]**2 * m.C.eta_o(f)) +
        m.C.Cr*ni**2 / (2 * m.C.M[0]**2 * m.C.eta_r(f))
        for i, f in enumerate(frange)
    ]
    for nv, ni in [
        [iNv(m.n0), iNi(m.n0)],
        [iNv(2*m.n0), iNi(2*m.n0)],
        [iNv(3*m.n0), iNi(3*m.n0)],
        [iNv(3.8*m.n0), iNi(3.8*m.n0)]
    ]
])
# print(frange.shape, Ueff_list.shape)
# U_out = np.insert(Ueff_list.transpose(), 0, frange, axis=1)
# print(U_out.shape)
# np.savetxt(join(m.foldername, 'U_eff.dat'), U_out, fmt='%.6f')
#
#
#
# plt.plot(wr.nrange[:len(n_V)]/wr.n0, n_V)
# plt.plot(wr.nrange[:len(n_I)]/wr.n0, n_I)
# plt.plot(wr.nrange[:len(n_I)]/wr.n0, Ueff)
# plt.show()


from Wrapper2 import Model

__author__ = 'const'
import Models2
from matplotlib import pyplot as plt
import numpy as np
import eosWrap as eos


# plt.plot(m.nrange, [m.sym.Ebind(m.nrange), m2.sym.Ebind(m.nrange)])
frange = np.linspace(0, 1, 100)

# plt.plot(frange, list(map(m.sym.C.eta_r, frange)))
# plt.plot(frange, list(map(m2.sym.C.eta_r, frange)))
# plt.show()
#
# plt.plot(frange, list(map(m.sym.C.eta_o, frange)))
# plt.plot(frange, list(map(m2.sym.C.eta_o, frange)))
# plt.show()
#
# plt.plot(frange, list(map(m.sym.C.eta_s, frange)))
# plt.plot(frange, list(map(m2.sym.C.eta_s, frange)))
# plt.show()
# E1, f1 = m.sym.Ebind(m.nrange, ret_f=1)
# E2, f2 = m2.sym.Ebind(m.nrange, ret_f=1)
# plt.plot(m.nrange, E1)
# plt.plot(m.nrange, E2)
# plt.show()
#
# plt.plot(m.nrange, f1)
# plt.plot(m.nrange, f2)
# plt.show()


wr = Models2.MKVOR_d()
m1 = wr.delta_phi
# wr.setDeltaConst(np.array([1.25, 1.25, 1.25, 1.25]),
#                  np.array([1., 1., 1., 1.]),
#                  np.array([1., 1., 1., 1.]),
#                  '1')

# wr.setDeltaConst(np.array([1.15, 1.15, 1.15, 1.15]),
#                  np.array([.9, .9, .9, .9]),
#                  np.array([1., 1., 1., 1.]),
#                  '2')
n = 12
print([m1.C.X_s[i] for i in range(n)])

print([m1.C.X_o[i] for i in range(n)])

print([m1.C.X_r[i] for i in range(n)])

print([m1.C.X_p[i] for i in range(n)])

print([m1.C.Q[i] for i in range(n)])

print([m1.C.S[i] for i in range(n)])
print([m1.C.T[i] for i in range(n)])

print(m2.C.f0)

wr.delta_phi.reset()
raise SystemExit
m1 = wr.delta
m2 = wr.delta_phi
m3 = wr.delta_phi_sigma


print([m1.C.X_s[i] for i in range(n)])

print([m1.C.X_o[i] for i in range(n)])

print([m1.C.X_r[i] for i in range(n)])

print([m1.C.X_p[i] for i in range(n)])

print([m1.C.Q[i] for i in range(n)])

print([m1.C.S[i] for i in range(n)])
print([m1.C.T[i] for i in range(n)])

print(m2.C.f0)

# exit()



for m in [m2]:
    for x in [[1.0, 1.0], [1.25, 1.], [1.15, .9]]:
        wr.setDeltaConst(np.array([x[0] for i in range(4)]),
             np.array([x[1] for i in range(4)]),
             np.array([1., 1., 1., 1.]),
             's=%.2f o=%.2f'%(x[0], x[1]))
        m.dumpEos()
        # m.dumpMassesCrust(nmax=m.nrange[m.rho.shape[0]-1])
        # m.dumpScalings()
        # m.dumpVs()
        # m.dumpChi()
        # m.dumpEtap()
        # m.dumpMeff()
        m.dumpMu()
        # m.dumpMassesCrust(nmin=2.4*m.n0, nmax=2.55*m.n0, fname='mass_crust_detail.dat')
        exit()

n = 12
print([m1.C.X_s[i] for i in range(n)])
print([m2.C.X_s[i] for i in range(n)])

print([m1.C.X_o[i] for i in range(n)])
print([m2.C.X_o[i] for i in range(n)])

print([m1.C.X_r[i] for i in range(n)])
print([m2.C.X_r[i] for i in range(n)])

print([m1.C.X_p[i] for i in range(n)])
print([m2.C.X_p[i] for i in range(n)])
#
# lw = 5
# m1.reset()
# fline = plt.plot(m1.nrange[:m1.concentrations().shape[0]]/m1.n0, m1.rho[:,0], lw=lw)
# lines = fline + plt.plot(m1.nrange[:m1.concentrations().shape[0]]/m1.n0, m1.concentrations(), lw=lw)
# plt.legend(lines, ['f'] + m1.part_names, loc=0, ncol=2)
# m2.reset()
# plt.gca().set_color_cycle(None)
# plt.plot(m2.nrange[:m2.concentrations().shape[0]]/m1.n0, m2.rho[:,0], ls='--', lw=lw)
# plt.plot(m2.nrange[:m2.concentrations().shape[0]]/m2.n0, m2.concentrations(), ls='--', lw=lw)
# font='32'
# plt.xlabel(u'$n/n_0$', fontsize=font)
# plt.ylabel(u'$n_i/n$', fontsize=font)
# plt.show()
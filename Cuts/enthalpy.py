from scipy.interpolate.interpolate import interp1d

__author__ = 'const'
from TOVInvert.Invert import EnthalpySolver
import Models2
from matplotlib import pyplot as plt

n_star = 2.5

m = Models2.KVOR()
m2 = Models2.KVORcut02()

E, P, n = m.nucl.EPN()
E2, P2, n2 = m2.nucl.EPN()

# plt.plot(E2, P2)
# plt.show()

P = P/m.mpi4_2_mevfm3
P2 = P2/ m.mpi4_2_mevfm3

iPofN = interp1d(n, P)
iPofN2 = interp1d(n2, P2)

out = m.nucl.dumpMasses(write=0, nmin=n_star, npoints=1, ret_frac=0)
out2 = m2.nucl.dumpMasses(write=0, nmin=n_star, npoints=1, ret_frac=0)

e_solver = EnthalpySolver(E[1:], P[1:])
e_solver2 = EnthalpySolver(E2[1:], P2[1:])

hc = .4

M, R, ph, eh, hrange = e_solver.integrateOut(e_solver.PofH(hc))
M2, R2, ph2, eh2, hrange2 = e_solver2.integrateOut(e_solver2.PofH(hc))

print(out[1], M[-1])
print(out[2], R[-1])

print(out2[1], M2[-1])
print(out2[2], R2[-1])

plt.plot(hrange, M)
plt.plot(hrange2, M2)
plt.show()

plt.plot(hrange, R)
plt.plot(hrange2, R2)
plt.show()

plt.plot(hrange, (M2 - M)/M, label='M')
plt.plot(hrange, (R2 - R)/R, label='R')
plt.legend()
plt.show()

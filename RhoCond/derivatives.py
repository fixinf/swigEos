import eosWrap as eos
import numpy as np
from matplotlib import pyplot as plt
from scipy.misc.common import derivative
import Models2

wr = Models2.KVOR()

m = wr.rcond_nucl
m2 = wr.nucl
m2.loadEos()
m.reset()
# m.dumpEos()
mus = m.mu()

dE_an = []

def func(f, n_n, n_p, mu_c):
	return eos.E_rho(np.array([f, n_n, n_p]), mu_c, m.C)

for i, n in enumerate(m.rho):
	ntot = n[1] + n[2]
	dEn = derivative(lambda z: func(n[0], z, n[2], m.mu_e[i]), n[1], dx=1e-3)
	dEp = derivative(lambda z: func(n[0], n[1], z, m.mu_e[i]), n[2], dx=1e-3)
	dE = n[1] / ntot * dEn + n[2]/ntot * dEp
	dE_an.append(dE)

plt.plot(mus[:, 0] - mus[:, 1])
plt.plot(m.mu_e)
plt.show()
# exit()

conc = m.concentrations()
l_conc = m.lepton_concentrations()
dE_mu = mus[:, 0] * conc[:, 0] + mus[:, 1] * conc[:, 1]
# plt.plot(m.nrange/m.n0, m.Efull())
plt.plot(m.nrange/m.n0, np.gradient(m.Efull(), m.nrange[1]-m.nrange[0]), label='grad')
plt.plot(m2.nrange/m2.n0, np.gradient(m2._E, m2.nrange[1]-m2.nrange[0]), label='old')
plt.plot(m.nrange/m.n0, dE_an, label='an')
plt.plot(m.nrange/m.n0, dE_mu, label='mu')
plt.legend()
plt.show()

# exit()

lines = plt.plot(m.nrange/m.n0, m.concentrations())
plt.plot(m2.nrange/m2.n0, m2.concentrations())
lines_rho = plt.plot(m.nrange/m.n0, m.nc/m.nrange)
plt.legend(lines+lines_rho, m.part_names + ['rho'])
plt.show()

plt.plot(m.nrange/m.n0, conc[:,1])
plt.plot(m.nrange/m.n0, l_conc[:, 0] + l_conc[:, 1] + m.nc/m.nrange)
plt.show()

plt.plot(m.nrange/m.n0, m.mu_e)
plt.plot(m.nrange/m.n0, np.gradient(m.mu_e, m.nrange[1]-m.nrange[0]))
plt.show()

# plt.plot(m.nrange/m.n0, mus)
# plt.plot(m.nrange/m.n0, np.gradient(mus[:, 0], m.nrange[1]-m.nrange[0]))
# plt.plot(m.nrange/m.n0, np.gradient(mus[:, 1], m.nrange[1]-m.nrange[0]))
# plt.show()

print(l_conc, l_conc.shape)

plt.plot(m.nrange, m.nrange)
plt.plot(m.nrange, m.rho[:, 1] + m.rho[:, 2])
plt.show()

sum1 = mus[:, 0] * m.rho[:, 1] + mus[:, 1] * m.rho[:, 2] + m.mu_e * m.nrange* l_conc[:, 0]+ m.nrange * m.mu_e*l_conc[:, 1]
sum2 = mus[:, 0] * m.nrange - m.mu_e * m.nc

plt.plot(m.nrange/m.n0, mus[:, 1] + m.mu_e)
plt.plot(m.nrange/m.n0, mus[:, 0])
plt.show()


plt.plot(m.nrange/m.n0, sum1, label='sum1')
plt.plot(m.nrange/m.n0, sum2, label='sum2')
plt.show()

pdiff = m.nrange*np.gradient(m.Efull(), m.nrange[1]-m.nrange[0]) - m.Efull()

plt.plot(m.nrange/m.n0, m.P_chem(), label='chem')

p_mu = m.nrange * mus[:, 0] - m.Efull()

plt.plot(m.nrange/m.n0, p_mu*m.mpi4_2_mevfm3, label='p_mu')
plt.plot(m.nrange/m.n0, m.P_chem2(), label='chem2')
plt.plot(m2.nrange/m2.n0, m2._P, label='orig')
plt.plot(m.nrange/m.n0, pdiff* m.mpi4_2_mevfm3, label='diff')
plt.legend()
plt.show()
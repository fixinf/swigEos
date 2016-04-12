import eosWrap as eos
import numpy as np
from numpy import pi, sqrt
from scipy.misc import derivative
from scipy import optimize

import Models2
from RhoCond import rho_wrap as rw
# wr = Models2.myMod()

# wr.hyper_phi_sigma.dumpEos()
# wr.hyper_phi_sigma.dumpMassesCrust()

wr = Models2.MKVOR2_fom_a(0.74, 0.0)

wr.delta_sym.dumpEos()
exit()
# # wr.hyper_phi.dumpEos()
# # wr.hyper_phi.dumpMassesCrust()

# # wr.hyper_phi_sigma.dumpEos()
# # wr.hyper_phi_sigma.dumpMassesCrust()
# wr.delta_phi.dumpEos()
# wr.delta_phi.dumpMassesCrust()

wr.delta_phi_sigma.dumpEos()
wr.delta_phi_sigma.dumpMassesCrust()

exit()

for fom in np.linspace(0.7, 0.9, 10):
	wr = Models2.MKVOR2_fom(fom)
	wr.dumpEos()
	wr.dumpProps()
	# exit()
	wr.dumpScalings()
	wr.delta_sym.dumpEos()
	print(wr.delta_sym.foldername, wr.delta_sym.filenames['eos'])
exit()
C = wr.hyper_phi.C

n = 5.
n_p = .5
mu_c = .5
f = 0.2
dx=1e-4

n_in = np.array([f, n-n_p, n_p, n_p, n_p, n_p, n_p, n_p])

for mu_c in np.linspace(0, 10, 100):
	print ('deriv - anal = ', eos.mu_deriv(n_in, 1, mu_c, C) - eos.mu_rho(n_in, 1, mu_c, C))
	# print (eos.mu(n_in, 1, C))
	# print ('anal = ', eos.mu_rho(n_in, 1, mu_c, C))

mu_n = derivative(lambda z: rw.E_r(z, n_p, f, mu_c, C), n-n_p, dx=dx)

print(mu_n)

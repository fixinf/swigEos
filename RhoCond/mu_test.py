import eosWrap as eos
import numpy as np
from numpy import pi, sqrt
from scipy.misc import derivative
from scipy import optimize

import Models2
from RhoCond import rho_wrap as rw

wr = Models2.KVOR()
C = wr.C

n = 5.
n_p = .5
mu_c = .5
f = 0.2
dx=1e-4

n_in = np.array([f, n-n_p, n_p])

print ('deriv = ', eos.mu_deriv(n_in, 1, mu_c, C))
# print (eos.mu(n_in, 1, C))
print ('anal = ', eos.mu_rho(n_in, 1, mu_c, C))

mu_n = derivative(lambda z: rw.E_r(z, n_p, f, mu_c, C), n-n_p, dx=dx)

print(mu_n)

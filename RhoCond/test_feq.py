import eosWrap as eos
from matplotlib import pyplot as plt
import numpy as np
from scipy.misc.common import derivative
import Models2


def eq_f(x, n, mu_c, C):
    def func(z):
        n_in = np.insert(n, 0, z)
        return eos.E_rho(n_in, mu_c, C)

    return [derivative(func, x, dx=1e-4),
            derivative(func, x, dx=1e-4, n=2)]

def eq_f_cpp(x, n, mu_c, C):
    return eos.wrap_func_feq_rho(x, n, mu_c, C)

n = np.array([2., .5])
mu_c = 0.

wr = Models2.KVOR()
C = wr.C

frange = np.linspace(0, 1, 500)

for mu_c in np.linspace(0, 1., 5):
    plt.plot(frange, [eq_f(f, n, mu_c, C) for f in frange])

plt.show()


for mu_c in np.linspace(0, 1., 5):
    plt.plot(frange, [eq_f_cpp(f, n, mu_c, C) for f in frange])

plt.show()


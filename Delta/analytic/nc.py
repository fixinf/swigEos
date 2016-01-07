from math import pi
from os.path import join
import matplotlib
from scipy import optimize
from scipy.interpolate.interpolate import interp1d
import numpy as np
from matplotlib import pyplot as plt
import Models2



matplotlib.rcParams.update({'font.size': 22})
matplotlib.rcParams.update({'lines.linewidth': 3})



m = Models2.Wal_d().delta_sym
C = m.C
fname_eos = '/home/const/MKVdelta_local/data_new/Wal_ds=1.31 o=1.00/eos_sym.dat'
fname_nc = '/home/const/MKVdelta_local/data_new/Wal_ds=1.31 o=1.00/nc.dat'

# fname_nc = '/home/const/MKVdelta_local/data_new/MKVOR_d/delta/DeltaSym/nc.dat'
# fname_eos = '/home/const/MKVdelta_local/data_new/MKVOR_d/delta/DeltaSym/eos_sym.dat'

# fname_nc = '/home/const/MKVdelta_local/data_new/MKVOR_ds=1.31 o=1.00/nc.dat'
# fname_eos = '/home/const/MKVdelta_local/data_new/MKVOR_ds=1.31 o=1.00/eos_sym.dat'
data = np.loadtxt(fname_eos, skiprows=1)
nc = np.loadtxt(fname_nc, skiprows=1)
nrange = data[:,0]
Yd_num = data[:, -1]
i_f = interp1d(data[:, 0], data[:, 1], bounds_error=0, fill_value=1.)
an_nc = []

for xs, n in nc:
    f = i_f(n)
    print(n, f)
    mn = C.M[0] * (1-f)
    md = C.M[10] * (1 - C.M[0]/C.M[10] * xs * f)
    _nc = 2 * (2 * mn * (md - mn)) ** (1.5) / (3*pi**2)
    an_nc.append(_nc/m.n0)

def nc_solve(xs):
    def eq(z, xs):
        f = i_f(z)
        mn = C.M[0] * (1-f)
        md = C.M[10] * (1 - C.M[0]/C.M[10] * xs * f)
        return (3 * pi**2 * z * m.n0 / 2) ** (2/3) / (2 * mn) - (md - mn)

    res = optimize.root(lambda z: eq(z, xs), 1.).x
    return res[0]

def test_nc():
    print(nc_solve(1.2))
    # exit()
    nc_sol = [nc_solve(z) for z in nc[:,0]]
    print(nc[:,1].tolist())
    print(an_nc)
    print(nc_sol)
    line1, =plt.plot(nc[:,0], nc[:, 1])
    line2, =plt.plot(nc[:,0], an_nc)
    line3, =plt.plot(nc[:,0], nc_sol)
    plt.legend([line1, line2, line3], ['exact', 'approx', 'approx_sol'])
    plt.show()

def Y_delta(n, xs):
    # nc = nc_solve(xs)
    nc = 2.31
    res = [0.]
    f = i_f(n)
    mn = C.M[0] * (1-f)
    md = C.M[10] * (1 - C.M[0]/C.M[10] * xs * f)
    if n > nc:
        res = optimize.root(lambda z: (n*m.n0 - z)**(2/3) - (mn/md)*(z/4)**(2/3) - (nc*m.n0)**(2/3), 0.).x
    return res[0]/m.n0/n

def Y_delta2(n, xs):
    nc = 2.31
    res = [0.]
    f = i_f(n)
    mn = C.M[0] * (1-f)
    md = C.M[10] * (1 - C.M[0]/C.M[10] * xs * f)
    if n > nc:
        res = [(4*(n**(2/3) - nc**(2/3))*(md/mn))**(3/2)]
    return res[0]/m.n0/n

def test_ndelta():
    Yd = np.array([Y_delta(n, 1.31) for n in nrange])
    Yd2 = np.array([Y_delta2(n, 1.31) for n in nrange])
    print(Yd)
    print(Yd2)
    plt.plot(nrange, Yd*nrange)
    plt.plot(nrange, Yd2*nrange)
    plt.plot(nrange, Yd_num*nrange)
    plt.plot(nrange, nrange-2.31)
    plt.show()

# test_nc()
test_ndelta()
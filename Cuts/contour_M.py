from scipy import optimize, interpolate
from scipy.cluster.vq import _krandinit
from scipy.interpolate.interpolate import interp1d
from Cuts.analytic_const import eq
import Models2
import numpy as np
# from analytic_const import eq
__author__ = 'const'

params = {
        'K0' : 250.,
        'f0' : 0.2,
        'Cs' : 196.3428132661,
        'Co' : 90.7681870142,
        'Cr' : 88.7261140316,
        'b' : 0.0089455122,
        'c' : 0.0077076578
    }
c = 0.2
m = Models2.Cubero_cut(c, K0=params['K0'], f0=params['f0'], J0=30., suffix='', params=params)
C = m.sym.C
mn = C.M[0]

def M(K,f0):
    res= eq([C.Cs/mn**2, C.b, C.c], [C.n0, -16./m.m_pi, mn, mn*(1-f0), K/m.m_pi, 30./m.m_pi])
    Cs = mn**2/res[0]
    b = res[1]
    c = res[2]
    Co = mn**2*res[3]
    Cr = mn**2 * res[4]
    m.setParams(Cs, Co, Cr, b, c, f0=f0)
    # print(m.sym.K())
    # print(m.nucl.Ebind(np.array([0, m.n0/2, m.n0])))
    m.nucl.reset()
    # n, mg, r, mb, mb2 = m.nucl.dumpMasses(nmin=.6, write=False, npoints=20, ret_frac=0)
    n, mg, r = m.nucl.dumpMasses(write=0, nmin=.6, nmax=5, npoints=20)
    i_m = interpolate.interp1d(n, -mg)
    res = optimize.minimize_scalar(i_m, bounds=[n[2], n[-2]], method='Bounded')
    # print(res)
    # print(-res.fun)
    return -res.fun

# M(400., 0.35)
#
# print(m.sym.K())
# print(m.sym.J())
# print(m.sym.Ebind(np.array([0, m.n0/2, m.n0]), ret_f=1))
# exit()


# M(250, 0.19)
frange = np.linspace(0.1, 0.35, 30)
Krange = np.linspace(200, 400, 30)
for c in [0.2]:
    m = Models2.Cubero_cut(c, K0=params['K0'], f0=params['f0'], J0=30., suffix='', params=params)
    # print(m.n0)
    # exit()
    C = m.sym.C
    table = []
    for f in frange:
        for K in Krange:
            print('K, f = ', [K, f])
            # M(K,f)
            # print('K = %.2f' %(m.sym.K()))
            mass = M(K, f)
            print('K = %.2f' %(m.sym.K()))
            n_pressure = np.linspace(2., 4.5, 6) * m.n0
            print(n_pressure/m.n0)
            # exit()
            m.sym.reset()
            es, ps, nrange = m.sym.EPN()
            plist = interp1d(nrange, ps)(n_pressure)
            table.append([1-f, K, mass] + list(plist))

    table = np.array(table)
    print(table)
    np.savetxt('contour_M_nocrust%.2f.dat'%(c), table)







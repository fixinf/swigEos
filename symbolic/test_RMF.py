from sympy.polys.groebnertools import is_reduced
from symbolic.RMF import RMF
import sympy as sp
import numpy as N
import Models2
from scipy.optimize import root
r = RMF()
M = Models2.waleckaMatsui()
print(r.Ek())
print(r.Em())
print(r.Etot())
print(r.Eqs()[1])
sp.var('go mo Co mn Cs nn np n L_o')
om = sp.Symbol('om')
om_sol = sp.solve(r.Eqs()[1], om)
om_real = sp.Mul()
for i in om_sol:
    if 'I' not in str(i):
        print('!')
        om_real = i
    print(i)
Etot = r.Etot().subs(om, om_real)
print(Etot.subs(go**2, Co / mn**2 * mo**2))
f = sp.Symbol('f')
print(r.Eqs()[0])
_eqfun = sp.lambdify([f, Cs, nn, np, mn], r.Eqs()[0], modules='numpy')

_Cs = 266.9
_Co = 195.7
_mn = M.C.M[0]
n0 = ((1.42*197.33/135)**3) / ((3*N.pi**2)/2.)

def f_eq(n):
    return root(lambda z: _eqfun(z, _Cs, n/2, n/2, _mn), 0.4).x[0]

print(f_eq(n0))
Etot = r.Etot().subs(om, om_real).subs(go, sp.sqrt(Co)/mn * mo)
print(Etot)
_Efun = sp.lambdify([f, Cs, Co, nn, np, mn, L_o], Etot, modules='numpy')
print((_Efun(f_eq(n0), _Cs, _Co, n0/2, n0/2, _mn, 0.01) / n0 - _mn)*135)

# _Efun = sp.lambdify([n])



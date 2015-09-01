__author__ = 'const'
import sympy as sp
import numpy as np
from sympy.utilities.lambdify import lambdify, lambdastr
import Models2
M = Models2.KVOR()
sp.var('Cs Co n x m f mn meff')
pf = (3 * sp.pi**2 * n)**(sp.Rational(1,3))
Ekin = 2 * (-m**4*sp.asinh(pf/m)/8 + m**3*pf/(8*sp.sqrt(1 + pf**2/m**2)) +
            3*m*pf**3/(8*sp.sqrt(1 + pf**2/m**2)) + pf**5/(4*m*sp.sqrt(1 + pf**2/m**2)))/sp.pi**2

E = mn**4 * f**2 / (2 * Cs) + Co * n**2 / (2 * mn**2) + Ekin.subs(m, mn*(1-f))

print(E)

_Cs = 266.9
_Co = 195.7
_mn = M.C.M[0]
n0 = ((1.42*197.33/135)**3) / ((3*np.pi**2)/2.)
# Estr = lambdastr([n, f, Cs, Co, meff, mn], E.doit())
# print(Estr)

Efun = lambdify([n, f, Cs, Co, meff, mn], E.doit(), modules='numpy')

Elambda = lambda n,f : Efun(n, f, _Cs, _Co, _mn * (1 - f), _mn)
print(135*(Elambda(n0, 0.44) / n0 - _mn))







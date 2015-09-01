__author__ = 'const'
import numpy as np
import sympy as sp

class RMF():
    def __init__(self):
        b_names = ['n', 'p']
        m_names = ['s', 'o', 'r']
        l_names = ['e', 'mu']
        self.mpi = 135.
        self.mn = 938/self.mpi
        b_masses = {'n' : self.mn, 'p': self.mn}
        pass

    def setupFunc(self):
        pass

    def Ek(self):
        sp.var('nn np n m mn f')
        pf = (3 * sp.pi**2 * n)**(sp.Rational(1, 3))
        int = (-m**4*sp.asinh(pf/m)/8 + m**3*pf/(8*sp.sqrt(1 + pf**2/m**2)) +
            3*m*pf**3/(8*sp.sqrt(1 + pf**2/m**2)) + pf**5/(4*m*sp.sqrt(1 + pf**2/m**2)))/sp.pi**2
        return [int.subs(m, mn * (1-f)).subs(n, nn), int.subs(m, mn * (1-f)).subs(n, np)]

    def Em(self):
        sp.var('nn f np mn Cs Co')
        return [mn**4 * f**2 / (2 * Cs), Co * (nn + np)**2 / (2 * mn**2)]

    def Etot(self):
        res = 0
        for _e in self.Ek():
            res = res + _e

        sp.var('mn go mo mn Cs f om L_o')
        E = res + mn**4 * f**2 / (2 * Cs) - mo**2 * om**2 / 2 + go * om * (nn + np) + L_o / 4 * om**4
        return E

    def Eqs(self):
        return [sp.diff(self.Etot(), f),
                sp.diff(self.Etot(), om)]
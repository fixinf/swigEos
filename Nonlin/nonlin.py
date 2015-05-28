__author__ = 'const'
import numpy as np
from scipy.optimize import root
from math import pi, asinh, sqrt


class NLVector():
    def __init__(self, Cs, Co, Cr, b, c, Co4, Cr4):
        self.Cs = Cs
        self.Co = Co
        self.Cr = Cr
        self.b = b
        self.c = c
        self.Co4 = Co4
        self.Cr4 = Cr4
        self.mpi = 135
        self.mn = 938/135
        self.mo = 783/135
        self.mr = 770/135
        self.go = sqrt(self.Co) /self.mn * self.mo
        self.gr = sqrt(self.Cr) /self.mn * self.mr
        self.n0 = ((1.42*197.33/135)**3) / ((3*np.pi**2)/2.)
        
    def pf(self, n):
        return (3 * pi**2 * n)**(1/3)

    def U(self, f):
        return 0.

    def I1(self, m, n):
        pf = self.pf(n)
        return (-m**4*asinh(pf/m)/8 + m**3*pf/(8*sqrt(1 + pf**2/m**2))
                + 3*m*pf**3/(8*sqrt(1 + pf**2/m**2)) +
                pf**5/(4*m*sqrt(1 + pf**2/m**2)))/pi**2


    def I2(self, m, n):
        pf = self.pf(n)
        return (-m**3*asinh(pf/m)/2 + m**2*pf/(2*sqrt(1 + pf**2/m**2)) +
                   3*pf**3/(8*sqrt(1 + pf**2/m**2)) + pf**3/(8*(1 + pf**2/m**2)**(3/2)) -
                   pf**5/(4*m**2*sqrt(1 + pf**2/m**2)) + 3*pf**5/(8*m**2*(1 + pf**2/m**2)**(3/2)) +
                   pf**7/(4*m**4*(1 + pf**2/m**2)**(3/2)))/pi**2



    def f_eq(self, nn, np):
        meff = lambda x: self.mn * (1 - x)
        eq = lambda z: self.mn**4 * z / (self.Cs) - self.mn * self.I2(meff(z), nn) - self.mn * self.I2(meff(z), np)\
                       + self.U(z)
        return root(eq, 0.5).x[0]

    def vec_eq(self, nn, np):
        eq = lambda x: [-self.mo**2 * x[0] + self.go * (nn + np) + self.Co4 * x[0]**3 + self.Cr4 * x[1]**2 * x[0],
                            -self.mr**2 * x[1] + self.gr * (np - nn)/2 + self.Cr4 * x[1] * x[0]**2]
        return root(eq, [1,1]).x

    def _E(self, f, om, rho, nn, np):
        E = self.mn**4 * f**2 / (2 * self.Cs)
        E += -self.mo**2 * om**2 /2 + self.go * om * (nn + np) + self.Co4 * om**4 / 4
        E += -self.mr**2 * rho**2 /2 + self.gr * rho * (np - nn)/4 + self.Cr4 * om**2 * rho**2
        E += self.I1(self.mn*(1-f), nn)
        E += self.I1(self.mn*(1-f), np)
        return E

    def E(self, nn, np):
        f = self.f_eq(nn, np)
        om, rho = self.vec_eq(nn, np)
        return self._E(f, om, rho, nn, np)




m = NLVector(266.9, 195.7, 0, 0, 0, 0, 0)
print(m.f_eq(m.n0/2, m.n0/2))
print(m.vec_eq(m.n0/2, m.n0/2))

print(135*(m.E(m.n0/2, m.n0/2)/m.n0 - m.mn))

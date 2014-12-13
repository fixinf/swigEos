import eosWrap as eos
from pylab import *
from scipy.integrate import quad

def I1(n, m):
    p = eos.p_f(n)
    res = sqrt(m**2 + p**2)*(-3*m**2 * p/8. + p**3 / 4.)
    res += 3 * log(p/m + sqrt(1 + p**2/m**2)) / 8.
    res /= (pi**2)
    return res

def I1num(n, m):
    p = eos.p_f(n)
    res = quad(lambda z: z**4/(z**2 + m**2)**(3./2), 0, p)[0]
    res /= pi**2
    return res


print I1(1., 1.)
print I1num(1., 1.)

def I2(n, m):
    p = eos.p_f(n)
    res = p * sqrt(p**2 + m**2)*(m**2 + 2*p**2)
    res -= m**4 * log(p/m + sqrt(1 + p**2 / m**2))
    res /= 8* pi **2
    return res

def I2num(n, m):
    p = eos.p_f(n)
    res = quad(lambda z: z**2*sqrt(z**2 + m**2), 0, p)[0]
    res /= pi**2
    return res 

def I3(n, m):
    p = eos.p_f(n)
    res = p * sqrt(p**2 + m**2)
    res -= m**2 * log(p/m + sqrt(1 + p**2 / m**2))
    res /= 2 * pi **2
    return res

def I3num(n, m):
    p = eos.p_f(n)
    res = quad(lambda z: z**2/sqrt(z**2 + m**2), 0, p)[0]
    res /= pi**2
    return res

def I4num(n, m):
    p = eos.p_f(n)
    res = quad(lambda z: z**2/(z**2 + m**2)**(3./2), 0, p)[0]
    res /= pi**2
    return res
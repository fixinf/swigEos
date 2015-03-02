import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join
from time import sleep
from scipy.misc.common import derivative
from pylab import sqrt, exp



C = Models.myMod()
wr = Wrapper(C)

# TYPEL

fname='typel_DD.dat'
G_s = 10.685257
n0 = 0.148746 * wr.n0 / 0.16
a_s = 1.371545
b_s = 0.644063
c_s = 1.034552
d_s = 0.567627

G_o = 13.312280

a_o = 1.385567
b_o = 0.521724
c_o = 0.869983
d_o = 0.618991

a_r = 0.4987

# DD-ME2

# fname='DD_ME2.dat'
# G_s = 10.685257
# n0 = 0.152 * wr.n0 / 0.16
# a_s = 1.3881
# b_s = 1.0943
# c_s = 1.7057
# d_s = 0.4421
#  
# G_o = 13.312280
#  
# a_o = 1.3892
# b_o = 0.9240
# c_o = 1.4620
# d_o = 0.4775
#  
# a_r = 0.5647

# DD-ME1

# fname='DD_ME1.dat'
# G_s = 10.685257
# n0 = 0.152 * wr.n0 / 0.16
# 
# a_s = 1.3854
# b_s = 0.9781
# c_s = 1.5342
# d_s = 0.4661
# 
# G_o = 13.312280
# 
# a_o = 1.3879
# b_o = 0.8525
# c_o = 1.3566
# d_o = 0.4957
# 
# a_r = 0.5008


# TYPEL DD-F

# fname='typel_DD_F.dat'
# G_s = 11.024
# n0 = 0.1469 * wr.n0 / 0.16
# a_s = 1.4867
# b_s = 0.19560
# c_s = 0.42817
# d_s = 0.88233
# 
# G_o = 13.575
# 
# a_o = 1.5449
# b_o = 0.18381
# c_o = 0.43969
# d_o = 0.87070
# 
# a_r = 0.44793
 


# DD-ME2

def Cs(n):
    x = n/n0
    res = a_s * (1 + b_s * (x + d_s)**2)/(1 + c_s * (x + d_s)**2)
    return res#*G_s

def Co(n):
    x = n/n0
    res = a_o * (1 + b_o * (x + d_o)**2)/(1 + c_o * (x + d_o)**2)
    return res#*G_o

def Cr(n):
    x = n/n0
    return exp(- a_r * (x - 1))

n = np.linspace(0., 8*n0, 80, endpoint=False)
print Cs(n0), Co(n0)
print derivative(Cs, 0., dx=1e-3, n=2), derivative(Co, 0., dx=1e-3, n=2)
print derivative(Cs, n0, dx=1e-3, n=2), derivative(Co, n0, dx=1e-3, n=2)

scaling = wr.dumpScalings(None)

chi_s = []
chi_o = []

for i, _n in enumerate(scaling[:,0]):
    chi_s.append(C.phi_n(0, scaling[i, 1]) / sqrt(scaling[i, 2]))
    chi_o.append(C.phi_n(0, scaling[i, 1]) / sqrt(scaling[i, 3]))

typel_tab = np.array([n/n0, Cs(n), Co(n), Cr(n)]).transpose()
typel_table = tabulate(typel_tab,['n/n_0', 'f_s', 'f_o', 'f_r'], 
                       tablefmt='plain')

with open(fname, 'w') as f:
    f.write(typel_table)

# print scaling[:,0]
l_sc = plt.plot(scaling[:, 0], scaling[:, 2], scaling[:, 0], scaling[:, 3])
l_chi = plt.plot(scaling[:, 0], chi_s, scaling[:, 0], chi_o)
l_typ = plt.plot(n/n0, Cs(n), n/n0, Co(n))
plt.legend(l_sc + l_chi + l_typ, [r'$\eta_\sigma$',r'$\eta_\omega$',
                                  r'$\frac{\Phi(f)}{\sqrt{\eta_\sigma}}$',
                                  r'$\frac{\Phi(f)}{\sqrt{\eta_\omega}}$',
                                  r'$f_\sigma^{TYP}$',
                                  r'$f_\omega^{TYP}$'],
           loc=0)
plt.show()

import eosWrap as eos
from Wrapper import Wrapper
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
# x = linspace(0.0, 3.0, 100)
# l = plot(x, sin(x),x, sin(2*x),x, sin(3*x))
# colors = [c.properties()['color'] for c in l]
# print colors
# show()

C = eos.KVOR()
C.SetHyperConstants(2)
C.sprime = 0
C.Csp = 1.0
print C.Cs, C.Co, C.Cr, C.b, C.c

print eos.EBind(np.array([0.195, 0.25, 0.25]), C)

print C.X_s[3]



# print eos.f_eq(np.array([0.25, 0.25, 0.0]), np.array([0.0]), 1, C)
wr = Wrapper(C)
wr.reset(hyper=1)
rho = []
for r in wr.rho:
    rho.append(r/sum(r))
    
plot(wr.n/wr.n0, rho)
show()
print eos.stepE(0.5, array([0.5]), array([1e-5]), 1, 300, C)
#

# print eos.f_eq(np.array([0.25, 0.25]),np.array([0.0]), 1, C)

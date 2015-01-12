import matplotlib
from scipy.misc.common import derivative
from time import sleep
from math import exp
from scipy.interpolate.interpolate import interp1d
from matplotlib.widgets import Slider
matplotlib.use('QT4Agg')
import numpy as np
import matplotlib.pyplot as plt
import Models
from Wrapper import Wrapper
from scipy import optimize
import eosWrap as eos

C = eos.InterpolatedScalings()
C.debug = 0
C1 = Models.KVOR()
wr = Wrapper(C)
wr1 = Wrapper(C1)


npoints = 100
mix = 0.8
n_target = 100
start = 0
stop = npoints/5
power_0 = 0
power_1 = 2

a_s = -0.09

frange = np.linspace(-1e-2, 1+1e-3, npoints)
etar = np.array(map(C1.eta_r, frange))
etao = np.array(map(C1.eta_o, frange))
u = np.array(map(C1.U, frange))
etas = np.array(map(C1.eta_s, frange))
# etas = 2*C.Cs*np.array(map(C1.U, frange))*frange**2/C.M[0]**4 + 1
# etas[0] = 1.

# etar[5:15]-=0.05
# plt.plot(frange, etas)
s_stop = stop
etas[0:s_stop] = 1 + a_s*(1 - frange[0:s_stop]/0.195)**4 
# plt.plot(frange, etas)
# plt.show()

C.set_eta_o(frange, etao)
C.set_eta_r(frange, etar)
C.set_eta_s(frange, etas)
C.set_U(frange, u)


n = np.linspace(0, wr.n0, 100)
data = np.loadtxt('/home/const/GrabbedFigures/HebelerSchwenk2014/HebelerSchwenk_Out.dat',
                   skiprows=1)

# plt.plot(frange, map(C.eta_r, frange), frange, map(C1.eta_r, frange))
# plt.show()



iLower = interp1d(data[:,0]*wr.n0, mix*data[:,2] + (1-mix)*data[:,1],
                  kind='linear')

_n = np.linspace(0, wr.n0, n_target)
target = iLower(_n)

def func(x):
    etar[start:stop] = x
    C.set_eta_r(frange, etar)

    res = (135*(wr.Eneutr(_n)/_n - C1.M[0]) - target)**2
    res[0] = 0
    sqsum = np.sum(res)
    return sqsum

def callback(x):
    print x

x = 0
def solve():
    print func(etar[start:stop])
    print etar[start:stop]
    res = optimize.minimize(func, etar[start:stop],options={'disp':1},
                      callback=callback, method='BFGS')
    
    global x
    x = res.x
    plt.plot(n/wr1.n0, iLower(n), c='black')
    plt.plot(data[:,0], data[:, 1:],c='black',lw='3')
    plt.plot(n/wr1.n0, 135*(wr.Eneutr(n)/n - C1.M[0]),c='blue')
    plt.show()
    
    plt.plot(frange, map(C.eta_r, frange))
    plt.plot(frange, map(C1.eta_r, frange))
    plt.xlim([0,0.2])
    plt.ylim([0,5])
    plt.show()
    
    
solve()
print x
# x = np.array([ 0.11275636 , 0.53558749 , 0.56976262 ,
#                0.58749857 , 0.62771884 , 0.68508924,
#                0.74709851 , 0.82433091  ,0.89304006 , 0.98390119])
# etar[start:stop] = x
C.set_eta_r(frange, etar)

# fsmall = np.linspace(0., 0.195, 20)
fsmall = frange[start:stop]
etar_old = np.array(map(C1.eta_r, fsmall))

etar_new0 = map(C.eta_r, fsmall)
# 
# fig, ax = plt.subplots()
# plt.subplots_adjust(left=0.25, bottom=0.25)
# ax.plot(fsmall, map(C.eta_r, fsmall))
# ax.plot(fsmall, map(C1.eta_r, fsmall))
#  
# ax_a = plt.axes([0.25, 0.1, 0.65, 0.03])
# sl_a = Slider(ax_a, 'A', -20, 20, valinit=0)
#  
# ax_b = plt.axes([0.25, 0.15, 0.65, 0.03])
# sl_b = Slider(ax_b, 'B', -20, 20, valinit=0)
#  
# ax_c = plt.axes([0.25, 0.2, 0.65, 0.03])
# sl_c = Slider(ax_c, 'C', -20, 20, valinit=0)
#  
# l, = ax.plot(fsmall, etar_old)
# a = 0
# b = 0
# c = 0
# def update(val):
#     global a, b, c
#     a = sl_a.val
#     b = sl_b.val
#     c = sl_c.val
# #     print a, b
#     l.set_ydata(etar_old - 
#                  a*(fsmall/C.f0)**3*(1-fsmall/C.f0)**3 + 
#                  b*(fsmall/C.f0)**4*(1-fsmall/C.f0)**3 + 
#                  c*(fsmall/C.f0)**3*(1-fsmall/C.f0)**4)
# sl_a.on_changed(update)
# sl_b.on_changed(update)
# sl_c.on_changed(update)
#  
#  
# plt.show()

C.f0 = 0.195


def f(x):
    print x
    a = x[0]
    b = x[1]
    c = x[2]
#     print fsmall/C.f0
    return np.sum(((etar_new0 - (etar_old - 
                 a*(fsmall/C.f0)**power_0*(1-fsmall/C.f0)**power_1 + 
                 b*(fsmall/C.f0)**(power_0+1)*(1-fsmall/C.f0)**power_1 + 
                 c*(fsmall/C.f0)**power_0*(1-fsmall/C.f0)**(power_1+1))))**2)

res = optimize.minimize(f, [0.5, 05., -0.1,], callback=callback, tol=1e-8,
                        method='BFGS')
print res
a, b, c = res.x 

etar1 = etar
etar_new = (etar_old - 
             a*(fsmall/C.f0)**power_0*(1-fsmall/C.f0)**power_1 + 
             b*(fsmall/C.f0)**(power_0+1)*(1-fsmall/C.f0)**power_1 + 
             c*(fsmall/C.f0)**power_0*(1-fsmall/C.f0)**(power_1+1))
etar1[start:stop] = etar_new
frange_dense = np.linspace(0., wr.n0, 100)
plt.plot(frange_dense, map(C.eta_r, frange_dense))
C.set_eta_r(frange, etar1)
plt.plot(frange_dense, map(C.eta_r, frange_dense))
plt.xlim([0, 0.2])
plt.ylim([0, 1])
plt.show()

plt.plot(n/wr1.n0, iLower(n), c='black')
plt.plot(data[:,0], data[:, 1:],c='black',lw='1')
plt.plot(n/wr1.n0, 135*(wr.Eneutr(n)/n - C1.M[0]),c='blue')
plt.show()

wr.dumpVs()
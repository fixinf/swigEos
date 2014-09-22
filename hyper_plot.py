import eosWrap as eos
from Wrapper import Wrapper
from pylab import *
from numpy import array, linspace
from matplotlib.widgets import Slider, CheckButtons, RadioButtons
# x = linspace(0.0, 3.0, 100)
# l = plot(x, sin(x),x, sin(2*x),x, sin(3*x))
# colors = [c.properties()['color'] for c in l]
# print colors
# show()

C = eos.KVOR_mod()
C.SetHyperConstants(2)

# C.alpha = 0.85
# C.z = 0.65
# 
# C.omega_a = 1.67
# C.omega_f = 0.44
# 
# C.rho_a = 0.0
# 
# C.phi_a = -1.01
# C.phi_f = 0.33
# 
# C.d = -4.21
# 
# C.phi_gamma = 3.0
# C.phi_z = 3.5
# 
# C.sprime = 0
# C.Csp = 380.0
# 
# C.beta = 4.54
# C.gamma = 3.78
# 
# 
# 
# f0 = 0.26
C.sprime = 0

C.Csp = 200

f0 = 0.26



wr = Wrapper(C)
for f in linspace(C.f0, f0, 20):
    C.f0 = f
    wr.solve(f0=C.f0)

# C.SetHyperConstants(2)
C.sprime=1
sp = 1 + C.sprime
C.phi_meson = 1

# wr.reset(hyper=1, nmax=3.5, npoints=2000)

UpperY = [15.639, 60, 140, 210] 
UpperX = [2, 2.75, 3.5, 4.5]
LowerX= [2, 2.5, 3.5, 4, 4.5]
LowerY = [7.5, 17.49, 40.0, 60.0, 60.0]

fig, ax = subplots(1,2)
n_p = linspace(0.0, 3.5, 1000)
p = wr.Psymm(n_p)

wr.C.SetHyperConstants(2)

maxes = []
rho = []
iter = 0
styles = ['-', '--', '-.']
colors = ['b', 'r', 'g', 'y', 'cyan', 'm', 'coral', 'black']
mlines=[]
ener = []

flist = []
phi = 1
for hyper in [1]:
    for sprime in [0,1]:#range(hyper+1):
        for E in [-30.0]:#, -30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0]:
            if 1:
                C.sprime = sprime
                sp = 1 + sprime
                ener.append(E)
                C.phi_meson = phi
#                 C.set_xs(array([1.0, 1.0, E, 30.0, 30.0, 30.0, 250.0, 250.0]))
                wr.reset(hyper=hyper, nmax=3.5, npoints=400, iter=100)
                for i in wr.rho:
                    rho.append(i[sp:]/sum(i[sp:]))
                    flist.append(i[:sp])
                    
                    
                if (sprime == 0 and hyper == 1):
                    ax[0].set_color_cycle(colors)
                    lines = ax[0].plot(wr.n/wr.n0, rho, linestyle=styles[iter])
                    ax[0].legend(lines, ['n','p',r'$\Lambda$',r'$\Sigma^-$',r'$\Sigma^0$',
                                         r'$\Sigma^+$', r'$\Xi^-$', r'$\Xi^0$'], fontsize = 10)  
                else:
                    ax[0].set_color_cycle(colors)
                    ax[0].plot(wr.n/wr.n0, rho, linestyle=styles[iter])
                 
                ax[0].set_xlabel(r'$n/n_0$', fontsize=14)
                ax[0].set_ylabel(r'$n_i / n$', fontsize=14)
                
                wr.setDriver()
                E, P, N = wr.EPN()
                print wr.rho
                print wr.rho.shape
                print E
                print E.shape, N.shape
                _n, _M, _R = wr.stars(0.5, 3.5, 100)
                mline, = ax[1].plot(_n/wr.n0, _M, linestyle=styles[iter])
                mlines.append(mline)
                ax[1].set_xlabel(r'$n/n_0$', fontsize=14)
                ax[1].set_ylabel(r'$M$', fontsize=14)
                print max(_M)
                maxes.append(max(_M)) 
                rho = []
             
        iter += 1
             
print maxes
labels = []
for i, m in enumerate(maxes):
    print m
    print ener[i]
    labels.append(r'$M_{max}=%.2f$' % (m))
ax[1].legend(mlines, labels, loc=0, fontsize=10)



# print eos.f_eq(array([wr.n0*0.17/0.16/2, wr.n0*0.17/0.16/2]), C)
# print eos.EBind(array([0.15, wr.n0*0.17/0.16/2, wr.n0*0.17/0.16/2]), C)
# print eos.J(wr.n0*0.17/0.16, C)
# print eos.K(wr.n0*0.17/0.16, C)
show()
flist = array(flist)
print flist.shape
print n_p.shape
plot(wr.n/wr.n0, flist)
show()

mlist = []
masses = []
for f in flist:
    masses = []
    for i in range(8):
        print i
        marg = C.X_s[i]*C.M[0]/C.M[i]*f[0] #+ C.X_sp[i]*C.M[0]/C.M[i]*f[1]
        if C.sprime:
            marg += C.X_sp[i]*C.M[0]/C.M[i]*f[1]
        masses.append(C.phi_n(marg))
    mlist.append(array(masses))

mlist = array(mlist)
print mlist.shape
lines = plot(wr.n/wr.n0, mlist)
legend(lines, ['n','p',r'$\Lambda$',r'$\Sigma^-$',r'$\Sigma^0$',
                                         r'$\Sigma^+$', r'$\Xi^-$', r'$\Xi^0$'])
show()
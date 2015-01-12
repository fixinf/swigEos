import eosWrap as eos
from Wrapper import Wrapper
from pylab import *
from numpy import array, linspace
from matplotlib.widgets import Slider, CheckButtons, RadioButtons
import Models
# x = linspace(0.0, 3.0, 100)
# l = plot(x, sin(x),x, sin(2*x),x, sin(3*x))
# colors = [c.properties()['color'] for c in l]
# print colors
# show()

C = Models.myMod()
C.rho_kind = 2
C.rho_f = 0.3
C.rho_a = 10
C.rho_val = 3.5
C.SetHyperConstants(2)
frange = np.linspace(0., 1., 1000)
plt.plot(frange, map(C.eta_r, frange))
plt.show()
# C.omega_a = 300.
# C.omega_f = 0.3
C.omega_kind = 1
# C.SetHyperConstants(2)
# C.alpha = 0.85
# C.omega_a = 12.0
# C.rho_a = 0.0
# C.phi_f = 0.35
# C.phi_a = -3.5
# C.omega_f = 0.55
# C.d = -10.0
# C.phi_gamma = -0.8
# C.z = 0.65
# C.phi_z = 3.1
# 
# f0 = 0.26

# 
# C.alpha = 1.1
# C.omega_a = 5.0
# C.rho_a = 0.0
# C.phi_f = 0.35
# C.phi_a = -1.5
# C.omega_f = 0.55
# C.d = -5.0
# C.phi_gamma = 2.0
# C.gamma = 3.0
# C.beta = 20.8
# f0 = 0.24
f0 = 0.195
print C.X_s[7]

print C.X_s[7]
# print C.phi_f, C.omega_f, C.rho_f


wr = Wrapper(C)

UpperY = [15.639, 60, 140, 210] 
UpperX = [2, 2.75, 3.5, 4.5]
LowerX= [2, 2.5, 3.5, 4, 4.5]
LowerY = [7.5, 17.49, 40.0, 60.0, 60.0]

fig, ax = subplots(2, 2)
n_p = linspace(0.0, 4.0, 1000)
p = wr.Psymm(n_p)
lp, = ax[0][0].semilogy(n_p[:]/wr.n0, p)
ax[0][0].plot(UpperX, UpperY, LowerX, LowerY)
ax[0][0].set_xlabel(r'$n/n_0$', fontsize=14)
ax[0][0].set_ylabel(r'$P_{symm}$', fontsize=14)

wr.reset(hyper=0, npoints=2000, iter=100)
wr.setDriver()

_n, _M, _R, mg= wr.stars(0.5, 4.0, 100)
mline, = ax[0][1].plot(_n/wr.n0, _M)
legend([mline],[max(_M)])

wr.C.SetHyperConstants(2)
print C.X_s[7]
# pause(2)
maxes = []
rho = []
iter = 0
styles = ['-', '--', '-.']
colors = ['b', 'r', 'g', 'y', 'cyan', 'm', 'k', 'black']
mlines=[]
ener = []
C.sprime = 0
C.csp = 400
for hyper in [0]:
    for phi in [0]:
        for E in [-30.0]:#, -30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0]:
            if 1:
                ener.append(E)
                C.phi_meson = phi
#                 C.set_xs(array([1.0, 1.0, E, 30.0, 30.0, 30.0, 250.0, 250.0]))
                wr.reset(hyper=hyper, npoints=1000, iter=100)
                for i in wr.rho:
                    rho.append(i[1:]/sum(i[1:]))
             
                if (phi == 0 and hyper == 1):
                    ax[1][0].set_color_cycle(colors)
                    lines = ax[1][0].plot(wr.n/wr.n0, rho, linestyle=styles[iter])
                    ax[1][0].legend(lines, ['n','p','L','S-','S0','S+', 'X-', 'X0'], fontsize = 12)  
                else:
                    ax[1][0].set_color_cycle(colors)
                    ax[1][0].plot(wr.n/wr.n0, rho, linestyle=styles[iter])
                    ax[1][0].plot(wr.n/wr.n0, [0.14 for i in wr.n], linestyle='--')
                 
                ax[1][0].set_xlabel(r'$n/n_0$', fontsize=14)
                ax[1][0].set_ylabel(r'$n_i / n$', fontsize=14)
                 
                wr.setDriver()
                E, P, N = wr.EPN()
                frange = np.linspace(0, 1., 100)
                print E.shape, N.shape
                ax[1][1].plot(frange, map(C.eta_r, frange), linestyle=styles[iter])
                ax[1][1].set_xlabel(r'$f$', fontsize=14)
                ax[1][1].set_ylabel(r'$\eta_\rho$', fontsize=14)
                _n, _M, _R, mg = wr.stars(0.5, 4.0, 100)
                mline, = ax[0][1].plot(_n/wr.n0, _M, linestyle=styles[iter])
                mlines.append(mline)
                ax[0][1].set_xlabel(r'$n/n_0$', fontsize=14)
                ax[0][1].set_ylabel(r'$M$', fontsize=14)
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
ax[0][1].legend(mlines, labels, loc=0, fontsize=10)



# print eos.f_eq(array([wr.n0*0.17/0.16/2, wr.n0*0.17/0.16/2]), C)
# print eos.EBind(array([0.15, wr.n0*0.17/0.16/2, wr.n0*0.17/0.16/2]), C)
# print eos.J(wr.n0*0.17/0.16, C)
# print eos.K(wr.n0*0.17/0.16, C)
show()
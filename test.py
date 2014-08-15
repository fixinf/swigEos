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

C.alpha = 0.85
C.omega_a = 12.0
C.rho_a = 0.0
C.phi_f = 0.35
C.phi_a = -3.5
C.omega_f = 0.55
C.d = -10.0
C.phi_gamma = -0.8
C.z = 0.65
C.phi_z = 3.1

f0 = 0.26

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
print C.X_s[7]

print C.X_s[7]
# print C.phi_f, C.omega_f, C.rho_f


wr = Wrapper(C)

for f in linspace(C.f0, f0, 200):
    C.f0 = f
    wr.solve(f0=C.f0) 

# pause(5)

UpperY = [15.639, 60, 140, 210] 
UpperX = [2, 2.75, 3.5, 4.5]
LowerX= [2, 2.5, 3.5, 4, 4.5]
LowerY = [7.5, 17.49, 40.0, 60.0, 60.0]

fig, ax = subplots(2,2)
n_p = linspace(0.0, 4.0, 1000)
p = wr.Psymm(n_p)
lp, = ax[0][0].semilogy(n_p[1:]/wr.n0, p)
ax[0][0].plot(UpperX, UpperY, LowerX, LowerY)
ax[0][0].set_xlabel(r'$n/n_0$', fontsize=14)
ax[0][0].set_ylabel(r'$P_{symm}$', fontsize=14)

wr.reset(hyper=0, npoints=1000, iter=100)
wr.setDriver()

_n, _M, _R = wr.stars(0.5, 4.0, 100)
mline, = ax[0][1].plot(_n/wr.n0, _M)


plt.subplots_adjust(bottom=0.35)

###Define Sliders

axOmA = plt.axes([0.1, 0.05, 0.2, 0.05])
slOmA = Slider(axOmA, r'$a_\omega$',-10.0, 20.0, valinit=C.omega_a)

axPhiA = plt.axes([0.1, 0.15, 0.2, 0.05])
slPhiA = Slider(axPhiA, r'$a_\Phi$', -10.0, 20.0, valinit=C.phi_a)

axD = plt.axes([0.1, 0.25, 0.2, 0.05])
slD = Slider(axD, r'$d$', -15.0, 15.0, valinit=C.d)

#Define checkbutton
axChM = plt.axes([0.35, 0.05, 0.2, 0.05])
chM = RadioButtons(axChM, 'M')

#Update function
needM = 1
def Update(val):
    C.omega_a = slOmA.val
    C.phi_a = slPhiA.val
    C.d = slD.val
    wr.solve()
    lp.set_ydata(wr.Psymm(n_p))
    if needM:
        wr.reset(hyper=0, npoints=1000, iter=100)
        wr.setDriver()
        _n, _M, _R = wr.stars(0.5, 4.0, 100)
        mline.set_ydata(_M)
    fig.canvas.draw_idle()

def Upd_M(val):
    needM = not needM

#Set Update on sliders' change

slOmA.on_changed(Update)
slPhiA.on_changed(Update)
# chM.on_clicked(Upd_M)
slD.on_changed(Update)

wr.C.SetHyperConstants(2)
C.set_xs(array([1.0, 1.0, -20.0, 30.0, 30.0, 30.0, 250.0, 250.0]))
print C.X_s[7]
# pause(2)
maxes = []
rho = []
iter = 0
styles = ['-', '--', '-.']
colors = ['b', 'r', 'g', 'y', 'cyan', 'm', 'k', 'black']
mlines=[]
ener = []
for hyper in []:
    for phi in range(hyper+1):
        for E in [-30.0]:#, -30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0]:
            if 1:
                ener.append(E)
                C.phi_meson = phi
                C.set_xs(array([1.0, 1.0, E, 30.0, 30.0, 30.0, 250.0, 250.0]))
                wr.reset(hyper=hyper, npoints=1000, iter=100)
                for i in wr.rho:
                    rho.append(i[1:]/sum(i[1:]))
            
                if (phi == 0 and hyper == 1):
                    ax[1][0].set_color_cycle(colors)
                    lines = ax[1][0].plot(wr.n/wr.n0, rho, linestyle=styles[iter], label=E)
                    ax[1][0].legend(lines, ['n','p','L','S-','S0','S+', 'X-', 'X0'], fontsize = 12)  
                else:
                    ax[1][0].set_color_cycle(colors)
                    ax[1][0].plot(wr.n/wr.n0, rho, linestyle=styles[iter])
                
                ax[1][0].set_xlabel(r'$n/n_0$', fontsize=14)
                ax[1][0].set_ylabel(r'$n_i / n$', fontsize=14)
                
                wr.setDriver()
                E, P, N = wr.EPN()
                print E.shape, N.shape
                ax[1][1].plot(N/wr.n0, E, linestyle=styles[iter])
                ax[1][1].set_xlabel(r'$n/n_0$', fontsize=14)
                ax[1][1].set_ylabel(r'$\varepsilon \quad [m_\pi^4]$', fontsize=14)
                _n, _M, _R = wr.stars(0.5, 4.0, 100)
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
    labels.append(r'$U_\Lambda=%.2f, M_{max}=%.2f$' % (ener[i], m))
ax[0][1].legend(mlines, labels, loc=0, fontsize=10)
ax[1][0].text(3, 0.4, r'$U_\Lambda=-40$')
ax[1][0].text(7, 0.2, r'$U_\Lambda=0$')


# print eos.f_eq(array([wr.n0*0.17/0.16/2, wr.n0*0.17/0.16/2]), C)
# print eos.EBind(array([0.15, wr.n0*0.17/0.16/2, wr.n0*0.17/0.16/2]), C)
# print eos.J(wr.n0*0.17/0.16, C)
# print eos.K(wr.n0*0.17/0.16, C)
show()
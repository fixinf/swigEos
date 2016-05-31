import eosWrap as eos
import Models2
import numpy as np
from matplotlib import pyplot as plt
from scipy.misc.common import derivative
from matplotlib.widgets import Slider
import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 3

def getParts(n, m, mu_e):
    f = np.linspace(0, 1, 100)
    C = m.C
    mn = C.M[0]
    dx = 1e-5
    order = 3
    part_f = C.M[0]**4 / (2 * C.Cs) * np.array([derivative(lambda z: z**2 * C.eta_s(z), _f, dx=dx, order=order) for _f in f]) 
    part_f += np.array([derivative(lambda z: C.U(z), _f, dx=dx) for _f in f])
    part_kin = []
    for i, _n in enumerate(n):
        part_kin.append([derivative(lambda z: eos.kineticInt(_n, C.M[i] * C.phi_n(i, z), 2.), _f, dx=dx, order=order) for _f in f])
    part_kin = np.array(part_kin)
    sum_kin = np.sum(part_kin, axis=0)

    n_V = np.sum([C.X_o[i] * _n for i, _n in enumerate(n)])
    res_om =  - np.array([derivative(lambda z: C.eta_o(z), _f, dx=dx, order=order) for _f in f]) * np.array([1.0/C.eta_o(_f)**2 for _f in f]) * n_V**2 * C.Co / (2 * C.M[0]**2)
    # plt.plot(f, np.array([derivative(lambda z: 1.0/C.eta_r(z), _f, dx=dx) for _f in f]))
    # plt.plot(f, np.array([1.0/C.eta_r(_f) for _f in f]))
    # plt.show()
    n_I = np.sum([C.X_r[i] * C.T[i] * _n for i, _n in enumerate(n)])
    res_rho = -np.array([derivative(lambda z: C.eta_r(z), _f, dx=dx, order=order) for _f in f]) *np.array([1.0/C.eta_r(_f)**2 for _f in f]) * n_I**2 * C.Cr / (2 * C.M[0]**2)
    
    parts_rcond = []

    mu_e = np.float64(mu_e)

    m_rho = np.float64(C.m_rho)

    n_rho = (2 * mn**2 * m_rho / C.Cr)  * np.array([
            C.eta_r(_f)**0.5 * C.phi_n(0, _f)**2 / np.float64(C.chi_prime(_f)) * (
                1 - mu_e / (m_rho * np.float64(C.phi_n(0, _f))) ) for _f in f
        ])

    eta_rp = lambda z: C.phi_n(0, z)**2 / C.chi_prime(z)**2

    n_rho_f = lambda z: (2 * mn**2 * m_rho / 
        C.Cr) * C.eta_r(z)**0.5 * C.phi_n(0, z)**2 / C.chi_prime(z) * (
                1 - mu_e / (m_rho * C.phi_n(0, z)))

    pr_1 = C.Cr * abs(n_I)/ (2 * mn**2) * np.array([derivative(
        lambda z: n_rho_f(z) / C.eta_r(z), _f, dx=dx, order=order) for _f in f])

    pr_2 = -C.Cr / (8 * mn**2 ) * np.array([derivative(
        lambda z: n_rho_f(z)**2 / C.eta_r(z), _f, dx=dx, order=order) for _f in f])

    for i, _n in enumerate(n_rho):
        if n_rho[i]/2 > abs(n_I):
            pr_1[i] = 0.
            pr_2[i] = 0.

    for i, _f in enumerate(f):
        if n_rho[i]/2 < abs(n_I):
            parts_rcond.append([
                m_rho * abs(n_I) * (1 - mu_e / (m_rho * C.phi_n(0, _f))) * (
                    - eta_rp(_f)**.5 / C.eta_r(_f)**.5), 

                m_rho * abs(n_I) * (1 - mu_e / (m_rho * C.phi_n(0, _f))) * (
                    C.phi_n(0, _f) * derivative(lambda z: eta_rp(z)**.5 / C.eta_r(z)**.5, _f, dx=dx, order=order)),

                - m_rho * abs(n_I) * mu_e / (m_rho * C.phi_n(0, _f)) * (eta_rp(_f)**.5 / C.eta_r(_f)**.5),

                - m_rho**2 * mn**2 * C.phi_n(0, _f)**2/ (2 * C.Cr) * (1 - mu_e / (m_rho * C.phi_n(0, _f)))**2 * derivative(
                    eta_rp, _f, dx=dx, order=order),

                2*m_rho**2 * mn**2 * eta_rp(_f) / (2 * C.Cr) * (1 - mu_e / (m_rho * C.phi_n(0, _f))) * (mu_e / (m_rho)),

                +m_rho**2 * mn**2 * eta_rp(_f) * C.phi_n(0, _f)/ ( C.Cr) * (1 - mu_e / (m_rho * C.phi_n(0, _f)))**2

                ])
            res_rho[i] = 0.
        else:
            parts_rcond.append([0., 0., 0., 0., 0., 0.])




    return f, part_f, part_kin, res_om, res_rho, parts_rcond, pr_1, pr_2, n_rho


wr = Models2.myMod()
m = wr.rcond_nucl
m.dumpEos()
f = m.rho[:, 0]
plt.plot(m.nrange/m.n0, m.mu_e)
plt.plot(m.nrange/m.n0, m.C.m_rho*(1 - f))
plt.show()
# m.inspect_f()
# m.loadEos()

n = [0.808321 * 7.48*wr.n0,  (1-0.808321)*7.48*wr.n0]
# n = np.array([4., 0])
# mu_e = 1.
mu_e = 232.452920 / m.m_pi

f, part_f, kin, om, rho, rcond, pr_1, pr_2, n_rho = getParts(n, m, mu_e)
rcond = np.array(rcond)

# plt.plot(f, part_f, label='f')
# for i, k in enumerate(kin):
#     plt.plot(f, kin[i], label='kin[%i]'%i)
# plt.plot(f, om, label='om')
# plt.plot(f, rho, label='rho')
# plt.legend()
# plt.ylim([-10, 10])
# plt.show()



# plt.plot(f, rcond[:, 0] + rcond[:, 1] + rcond[:, 2], label='sum')
# plt.plot(f, pr_1, label='pr1')
# plt.legend()
# plt.ylim([-10, 10])
# plt.show()


# plt.plot(f, rcond[:, 3] + rcond[:, 4] + rcond[:, 5], label='sum')
# plt.plot(f, pr_2, label='pr2')
# plt.legend()
# # plt.ylim([-20, 20])
# plt.show()

# print(rcond.shape)
# # exit()

# sum_rcond = np.sum(rcond, axis=1)

# plt.plot(f, part_f + kin[0] + kin[1] + om + rho + sum_rcond, label='parts')

# plt.plot(f, part_f + kin[0] + kin[1] + om + rho + pr_1 + pr_2, label='parts2')
# plt.plot(f, [m.func_f(z, m.swArr(n), mu_e) for z in f], label='func')
# plt.plot(f, pr_1, label='pr1')
# plt.plot(f, pr_2, label='pr2')
# for i, _r in enumerate(rcond.transpose()):
#     plt.plot(f, _r, label="%i"%i)
# plt.plot(f, rho, label='rho')
# plt.plot(f, [0. for _f in f])
# plt.ylim([-10, 10])
# plt.legend()
# plt.show()


fig, ax = plt.subplots()

lf, = ax.plot(f, [m.func_f(z, m.swArr(n), mu_e) for z in f], label='func')
# lf ,= plt.plot(f, part_f + kin[0] + kin[1] + om +, label='parts2')
lpr1, = ax.plot(f, pr_1, label='pr1')
lpr2, = ax.plot(f, pr_2, label='pr2')
lr,= ax.plot(f, rho, label='rho')
lo, = ax.plot(f, om, label='om')
lnrho, = ax.plot(f, n_rho, label='n_rho')

ax.plot(f, [0. for _f in f])
ax.legend()
plt.subplots_adjust(left=0.25, bottom=0.25)
axn = plt.axes([0.25, 0.1, 0.65, 0.03])
sn = Slider(axn, 'N', 0, 10., valinit=n[0] / wr.n0)

axmu = plt.axes([0.25, 0.15, 0.65, 0.03])
smu = Slider(axmu, 'mu', 0, 10., valinit=mu_e)


def update(val):
    n = np.array([sn.val * wr.n0, 0])
    mu_e = smu.val
    f, part_f, kin, om, rho, rcond, pr_1, pr_2, n_rho = getParts(n, m, mu_e)
    lf.set_ydata([m.func_f(z, m.swArr(n), mu_e) for z in f])
    lpr1.set_ydata(pr_1)
    lpr2.set_ydata(pr_2)
    lr.set_ydata(rho)
    lo.set_ydata(om)
    lnrho.set_ydata(n_rho)
    fig.canvas.draw_idle()

sn.on_changed(update)
smu.on_changed(update)
ax.set_ylim([-20, 20.])
plt.show()
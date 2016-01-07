import eosWrap as eos
from scipy.interpolate.interpolate import interp1d
import Models2
import numpy as np
from matplotlib import pyplot as plt



wr = Models2.KVOR()
wr.loadEos()
m = wr.nucl
m.reset()
m.setCrust(inter='cubic', ncut_crust=.45, ncut_eos=.7)
exit()
E, P, N = m.EPN()
## mu = mu_n * x_n + mu_p * x_p

mu = np.nan_to_num((P + E) / N)

e = []
p = []
n = []
ncut_crust = 5
with open("/home/const/workspace/swigEosWrapper/crust.dat", 'r') as f:
    for line in f:
        # print line
        _e, _p, _n = line.split()
        if float(_n) < ncut_crust:
            e.append(float(_e))
            p.append(float(_p))
            n.append(float(_n))
e = np.array(e)/m.m_pi**4
p = np.array(p)/m.m_pi**4
n = np.array(n)*m.n0
# print(e.shape, p.shape, n.shape)
mu_crust = np.nan_to_num( (p + e) / (n))


npoints = 5000
murange = np.linspace(min(mu[1], mu_crust[1]), mu_crust[-1], npoints)
print(mu, mu_crust, murange[0])
iP = interp1d(mu, P)
ip = interp1d(mu_crust, p)

Pmu = iP(murange)
pmu = ip(murange)

i_cross = np.argmin(abs(pmu[3:] - Pmu[3:]))

width = int(0.03*npoints)

# num = int(0.02*npoints)
num = 2
print (i_cross, width, num)

mu_in = np.concatenate((murange[i_cross - width - num:i_cross-width],
                                         murange[i_cross + width : i_cross + width + num]))

p_in = np.concatenate((pmu[i_cross - width - num : i_cross - width],
                                         Pmu[i_cross + width : i_cross + width + num]))

print(mu_in)
print(p_in)

inter_coeff = np.polyfit(mu_in, p_in, 2*num)

pnew = np.concatenate((pmu[0:i_cross-width-num],
                       np.poly1d(inter_coeff)(murange[i_cross-width-num:i_cross+width+num]),
                       Pmu[i_cross+width+num:]))
print(murange.shape, pnew.shape)
# plt.plot(murange, np.poly1d(inter_coeff)(murange))
lines = plt.plot(murange, pmu, murange, Pmu)
plt.plot(murange, pnew)
plt.legend(lines, ['crust', 'bar'])
# l_mu_b ,= plt.plot(mu, P)
# l_mu_c ,= plt.plot(mu_crust, p)
# plt.legend([l_mu_b, l_mu_c], ['bar', 'crust'])
plt.show()


lines = plt.plot(N/m.n0, P, n/m.n0, p)
plt.legend(lines, ['bar', 'crust'])
plt.show()

lines = plt.plot(mu, N, mu_crust, n)
plt.legend(lines, ['bar', 'crust'])
plt.show()

m.reset()
mu_n = m.mu()[:,0]

lines = plt.plot(N, mu, n, mu_crust, m.nrange, mu_n)
plt.legend(lines, ['bar', 'crust'])
plt.show()
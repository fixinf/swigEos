import eosWrap as eos
import Models2
import numpy as np

wr = Models2.MKVOR2_fom(0.6)
wr.rcond_nucl.loadEos()
wr.rcond_nucl.getFnBranches()
exit()



m1 = wr.rcond_hyper_phi
m2 = wr.rcond_delta_phi



i = 1000
print("ntot = %.3f" % (m1.nrange[i]/wr.n0))
n = np.insert(m1.rho[i, :], len(m1.rho[i, :]), [m1.mu_e[i]])
# n = np.insert(m1.rho[i, :], len(m1.rho[i, :]), [0., 0., 0., 0., m1.mu_e[i]])
print(n)

n2 = np.insert(m2.rho[i, :], len(m2.rho[i, :]), [m2.mu_e[i]])
print(n2)
print('\n')
# print(m2.mu_e)

res = eos.wrap_fun_rho(n, m1.C, len(n)-2)
res2 = eos.wrap_fun_rho(n2, m2.C, len(n2)-2)
print(res) 
print(res2)
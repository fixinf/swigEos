__author__ = 'const'
import eosWrap as eos
import Models2

dr = eos.KVDriver()
m = Models2.KVOR()
m.nucl.check()
E, P, n = m.nucl.EPN()
P /= m.mpi4_2_mevfm3
# E *= m.m_pi ** 4
# P *= m.m_pi ** 4
dr.set(E, P, n)

res = eos.star_crust_i(1.5, 4, dr, 1e-11)
print(res)
print(res[3]/1.4766 * 2e33 * 1e10)


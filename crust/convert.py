#!/usr/bin/python
from Wrapper import Wrapper
import eosWrap as eos

C = eos.KVOR()
wr = Wrapper(C)

eps = []
p = []
nb = []

with open("ex_bps_input.txt", 'r') as f:
	for i,line in enumerate(f):
		if i > 0:
			print line
			_eps, _p, _nb = line.split()
			eps.append(float(_eps))
			p.append(float(_p))
			nb.append(float(_nb))
print eps
print p
print nb

eps_out = []
p_out = []
nb_out = []

eps_mult = 4.31e-6 / wr.m_pi**4
p_mult = 4.795e-27 / wr.m_pi**4
nb_mult = 1e-39/0.16

with open("crust.dat", 'w') as f:
	for i, _i in enumerate(nb):
		eps_out.append(eps_mult*eps[i])
		p_out.append(p[i]*p_mult)
		nb_out.append(nb[i]*nb_mult)
		t_line = (eps_out[i], p_out[i], nb_out[i])
		print t_line
		f.write('%e %e %e \n' % t_line)






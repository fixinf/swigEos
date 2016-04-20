import Models2

wr = Models2.MKVOR2final()

wr.rcond_delta_phi.dumpEos()
wr.rcond_delta_phi.dumpMassesCrust()

wr.rcond_delta_phi_sigma.dumpEos()
wr.rcond_delta_phi_sigma.dumpMassesCrust()
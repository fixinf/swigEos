import Models2

wr = Models2.MKVOR2final()



for m in [wr.nucl, wr.hyper_phi, wr.hyper_phi_sigma,
            wr.delta_phi, wr.delta_phi_sigma]:
    m.loadEos()
    m.dumpScalingsN()

for m in [wr.rcond_nucl, wr.rcond_hyper_phi, wr.rcond_hyper_phi_sigma,
            wr.rcond_delta_phi, wr.rcond_delta_phi_sigma]:
    m.loadEos()
    m.dumpScalingsN()


for m in [wr.rcc_nucl, wr.rcc_hyper, wr.rcc_hyper_phi, wr.rcc_hyper_phi_sigma, 
            wr.rcc_delta_phi, wr.rcc_delta_phi_sigma]:
    m.loadEos()
    m.dumpScalingsN()
    # m.dumpMassesCrust()
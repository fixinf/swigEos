import Models2

wr = Models2.MKVOR_tanh2_cut(50., 0.6, 4)
wr.setDeltaPotential(-50.)

m = wr.rcond_delta_phi
m.loadEos()

if m.needsMaxw():
    m.processMaxw()

out = m.dumpMassesCrust()
print(out[1])

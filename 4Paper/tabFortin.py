__author__ = 'const'
import Models2

models = [Models2.KVOR(), Models2.KVORcut04(), Models2.KVORcut03(),
          Models2.KVORcut02(), Models2.myMod()]
# models = [Models2.KVOR()]


for m in models:
    # for s in [m.nucl, m.hyper, m.hyper_phi, m.hyper_phi_sigma]:
    for s in [m.nucl]:
        s.dumpPf()

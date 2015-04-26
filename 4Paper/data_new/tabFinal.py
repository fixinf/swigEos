__author__ = 'const'
import Models2
models = [Models2.KVOR(), Models2.myMod(), Models2.KVORcut02(),
          Models2.KVORcut03(), Models2.KVORcut04()]
for m in models:
    m.dumpEos()


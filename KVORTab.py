import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join
from pylab import pi, sin

model = Models.myMod
C = model()
C.set_xs(np.array([0., 0., -28., 30., 30., 30., -15., -15.]))
wr = Wrapper(C)
print C.X_s[2], C.X_s[3], C.X_s[6]
exit()
# wr.testHyperBind()
# J = wr.J()
# L = wr.L()
# Ksymm = wr.Ksymm()
# Kprime = wr.Kprime()
# K = wr.K() 
# print 'J=', J
# print 'L=', L
# print 'Ksymm=',Ksymm
# 
# print 'K=', K,
# print "K'=", Kprime
#   
# print 'eq(5):' , L, -11.76 + 3 * J + Ksymm/4.55
# print 'eq(8):', Ksymm, -307.862 + 3.292 * L
# exit()


# 
# nU = np.array([0.25, wr.n0/2])
# lines = []
# for i in range(8):
#     e, u = wr.UofE(i, nU)
#     line, = plt.plot(e, u)
#     lines.append(line)
# plt.legend(lines, ['n', 'p', 'L', 'Sm', 'S0', 'S+', 'X-', 'X0'])
# plt.show()
# exit()
folderName = join('/home/const/Dropbox/Documents/For DN/Very final/data',
                  model.__name__)
print folderName

wr.testPodsiedlowski(0.6, 0.8, folderName)

wr.dumpUofE(folderName)
exit()

# wr.dumpPotentialsNS(folderName, show=True)
# wr.dumpHinHPotentials()
# exit()
wr.dumpLandauParams(folderName)
exit()
wr.dumpHyper(folderName)
exit()
wr.dumpAll(folderName, folderName+'.zip')
# wr.testPodsiedlowski(0.8, 0.9, folderName)
# wr.dumpEos(folderName)
# wr.dumpMasses(folderName)
wr.dumpMassesCrust(ncut_crust=0.6, ncut_eos=0.8, folderName=folderName)
# wr.dumpHyper(folderName)




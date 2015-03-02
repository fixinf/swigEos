import matplotlib
from scipy.misc.common import derivative
from time import sleep
from math import exp
from cmath import sqrt
matplotlib.use('QT4Agg')
from matplotlib.widgets import Slider
import numpy as np
import matplotlib.pyplot as plt
import Models
from Wrapper import Wrapper
from scipy import optimize

C = Models.myModFinal()
C1 = Models.myMod()
wr = Wrapper(C)
wr1 = Wrapper(C1)

# wr.reset()
# wr1.reset()

wr.testDanielewicz()
wr1.testDanielewicz()

wr.reset()
wr1.reset()

plt.plot(wr.n/wr.n0, wr.concentrations(), wr1.n/wr1.n0,wr1.concentrations())
plt.show()
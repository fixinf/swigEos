'''
Created on 15 авг. 2014 г.

@author: const
'''
import eosWrap as eos
from Wrapper import Wrapper
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
# x = linspace(0.0, 3.0, 100)
# l = plot(x, sin(x),x, sin(2*x),x, sin(3*x))
# colors = [c.properties()['color'] for c in l]
# print colors
# show()

C = eos.KVOR()
C.SetHyperConstants(2)




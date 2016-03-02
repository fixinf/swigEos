from os.path import join

from scipy.interpolate.interpolate import interp1d

import Models2

__author__ = 'const'
import matplotlib
matplotlib.use("QT5Agg")
from matplotlib import  pyplot as plt
from RhoCond.rho_wrap import *
import numpy as np


wr = Models2.KVOR()


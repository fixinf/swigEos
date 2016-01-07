__author__ = 'const'
import numpy as np

mpi = 139.
mn = 938./mpi
md = 1232/mpi

n0 = .16 * (197.33/135.)**3

def Phi(n):
    return 1 - (n / 8 / n0)

def E(nn, nd):
    n = nn + nd
    return mn * nn * Phi(n) + md * nd * Phi(n) +
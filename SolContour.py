import numpy as np
from matplotlib import pyplot as plt
import eosWrap as eos

class SolContour:
    def __init__(self, m, axes_keys=['n', 'f'], **kwargs):
        for key, val in kwargs.items():
            setattr(self, key, val)
        if not set(axes_keys).issubset(set(kwargs.keys())):
            raise ValueError('Axes keys should be in the kwargs!')
        self.axes_keys = axes_keys
        self.model = m
    
    def plotArgsContour(self, var):
        val = getattr(self, var)
        X, Y = np.meshgrid(self.n, self.f)
        return X, Y, val

    def plotArgsLine(self, var):
        return self.n, getattr(self, var)
            
    def getLineConcs(self, x, y):
        ix = [np.argmin(abs(_x - getattr(self, self.axes_keys[0]))) for _x in x]
        iy = [np.argmin(abs(_y - getattr(self, self.axes_keys[1]))) for _y in y]    
        concs = np.array([self.rhos[_ix][_iy] for _ix, _iy in zip(ix, iy)])
        return concs

    def getLineMus(self, x, y):
        ix = [np.argmin(abs(_x - getattr(self, self.axes_keys[0]))) for _x in x]
        iy = [np.argmin(abs(_y - getattr(self, self.axes_keys[1]))) for _y in y]
        mues = np.array([self.mue[_ix][_iy] for _ix, _iy in zip(ix, iy)])
        return mues

    def __repr__(self):
        s = ''
        for k in self.__dict__.keys():
            try:
                s += '    |->%s: %s \n' % (k, getattr(self, k).shape.__repr__())
            except AttributeError:
                pass
        return "SolContour:\n" + s

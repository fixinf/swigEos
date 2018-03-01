import matplotlib as mpl
import Label
import matplotlib.pyplot as plt

standard_lss = ['-', '--', '-.', ':']

def ConcFnPlot(models, lss=standard_lss, num_legends=2):
    flag = 0
    lmod = []
    for m, ls in zip(models, lss):
        assert m.set
        if not flag:
            lcs = plt.plot(m.nrange/m.n0, m.concentrations())
            lf = plt.plot(m.nrange/m.n0, m.rho[:, 0], c='black')
            flag = 1
            _lcs = lcs
            _lf = lf
        else: 
            _lcs = [plt.plot(m.nrange/m.n0, conc, ls=ls, c=l.get_c())
            for l, conc in zip(lcs, m.concentrations().transpose())]
            _lf = plt.plot(m.nrange/m.n0, m.rho[:, 0], ls=ls, c=lf[-1].get_c())
        
        lmod += [_lcs[0]]

def ConcFnNcPlot(models, lss=standard_lss, num_legends=2):
    flag = 0
    lmod = []
    for m, ls in zip(models, lss):
        assert m.set
        if not flag:
            lcs = plt.plot(m.nrange/m.n0, m.concentrations())
            lf = plt.plot(m.nrange/m.n0, m.rho[:, 0], c='black')
            lnc = plt.plot(m.nrange/m.n0, m.nc/m.nrange, c='red')
            flag = 1
            _lcs = lcs
            _lf = lf
            _lnc = lnc
        else: 
            _lcs = [plt.plot(m.nrange/m.n0, conc, ls=ls, c=l.get_c())
            for l, conc in zip(lcs, m.concentrations().transpose())]
            _lf = plt.plot(m.nrange/m.n0, m.rho[:, 0], ls=ls, c=lf[-1].get_c())
            _lnc = plt.plot(m.nrange/m.n0, m.nc/m.nrange, c=lnc[-1].get_c())
        
        lmod += [_lcs[0]]
        return lmod
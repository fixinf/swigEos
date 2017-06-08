import matplotlib
from types import ModuleType

def get_axes(obj):
    # checks if the obj is pyplot
        if isinstance(obj, ModuleType):
            return obj.gca()
        else:
            return obj

def set_labels(obj, xlabel=None, ylabel=None):
    ax = get_axes(obj)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

def set_lim(obj, xlim=None, ylim=None):
    ax = get_axes(obj)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

def P_n(obj):
    ylabel = r'$P(n)$ [MeV/fm$^3$]'
    xlabel = r'$n/n_0$'
    set_labels(obj, xlabel, ylabel)

def D_f(obj):
    ylabel = r'$D(f), \, m_\pi^4$'
    xlabel = r'$f$'
    ylim = [-10, 10]
    xlim = [0, 1]
    set_labels(obj, xlabel, ylabel)
    set_lim(obj, xlim, ylim)
    
def Conc_n(obj):
    ylabel = r'$n_i/n$'
    xlabel = r'$n/n_0$'
    ylim = [0, 1]
    xlim = None
    set_labels(obj, xlabel, ylabel)
    set_lim(obj, xlim, ylim)

def M_n(obj):
    ylabel = r'$M/M_\odot$'
    xlabel = r'$n_{\rm cen}/n_0$'
    ylim = None
    xlim = None
    set_labels(obj, xlabel, ylabel)
    set_lim(obj, xlim, ylim)
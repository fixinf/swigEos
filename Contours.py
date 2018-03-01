import numpy as np
from matplotlib import pyplot as plt
from numpy import linspace
from SolContour import SolContour
import eosWrap as eos
import fixedNF
import Wrapper2
from Tqdm_wrap import tqdm

from scipy.misc.common import derivative
wr = Wrapper2.Model()

def getContour(m=wr, n_iter=300, 
               nrange=linspace(0, 8, 100), frange=linspace(0, 0.95, 100)):
    Ess = []
    Dss = []
    concss = []
    muess = []
    Rss = []
    for n in tqdm(nrange):
        i_init = np.argmin(abs(m.nrange/m.n0 - n))
        init = np.concatenate((m.rho[i_init, 2:], [m.mu_e[i_init]]))
        Es = []
        Ds = []
        mues = []
        concs = []
        Rs = []
        for f in frange:
            res, E = m.E_NF(n*wr.n0, f, init, iterations=n_iter)
            init = res[:-1].copy()
            n_n = n*wr.n0 - sum(init[:-1])
            conc = np.insert(init[:-1], 0, n_n)
            concs.append(conc)
            Es.append(E[0])
            mues.append(init[-1])
            D = eos.wrap_func_feq_rho(f, conc, init[-1], m.C)
#             print(D)
            Ds.append(D)
            Rs.append(res[-1])
        Ess.append(Es)
        concss.append(concs)
        Dss.append(Ds)
        muess.append(mues)
        Rss.append(Rs)
    return SolContour(m,
        f=frange,
        n=nrange,
        E=np.array(Ess).transpose(),
        D=np.array(Dss).transpose(),
        rhos=np.array(concss),
        mue=np.array(muess).transpose(),
        nc = np.array(Rss).transpose()
    )

def getZeroLine(cs, n_iter=100):
    plt.ioff()
    ct = plt.contour(*cs.plotArgsContour('D'), levels=[0.])
    paths = ct.collections[0].get_paths()
    plt.close()
    plt.ion()
    branches = []
    for p in paths:
        x, y = p.vertices.transpose()
        
        concs = cs.getLineConcs(x, y)
        mue = cs.getLineMus(x, y)
        
        newE = []
        newC = []
        newMu = []
        newD = []
        newNc = []
        newdD = []
        m = cs.model
        for n, f, mu, conc in zip(tqdm(x), y, mue, concs):
            init = np.concatenate((conc[1:], [mu]))
            res, E = m.E_NF(n*m.n0, f, init, iterations=n_iter)
            init = res[:-1].copy()
            n_n = n*m.n0 - sum(init[:-1])
            conc = np.insert(init[:-1], 0, n_n)
            newC.append(conc)
            newE.append(E[0])
            newMu.append(init[-1])
            newNc.append(res[-1])
            D = eos.wrap_func_feq_rho(f, conc, init[-1], m.C)
    #             print(D)
            newD.append(D)

            dD = derivative(lambda z: fixedNF.D_consistent(m, n, z, init, m.C, iterations=n_iter),
                f, dx=1e-3, order=3)
            newdD.append(dD)
            
        branches.append(SolContour(m,
            f = y,
            n = x,
            E=np.array(newE).transpose(),
            D=np.array(newD).transpose(),
            dD=np.array(dD).transpose(),
            rhos=np.array(newC),
            mue=np.array(newMu).transpose(),
            nc=np.array(newNc).transpose()
        ))
    return branches

class InterpMin:
    def __init__(self, E):
        self.inters = [
            interp1d(b.n, b.E) for b in branches
        ]
    
    def __call__(self, n):
        return min([i(n) for i in self.inters])




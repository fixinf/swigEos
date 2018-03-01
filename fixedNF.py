import numpy as np
import eosWrap as eos
from scipy.misc.common import derivative
from Tqdm_wrap import tqdm

def get_Dnf(m, nlist=np.linspace(0, 8, 100), 
            frange=np.linspace(0, 1, 100)):
    if not hasattr(m, 'rho'):
        raise ValueError('Reset the model explicitly to provide initial values!')

    Ess = []
    Dss = []
    Rss = []
    concss = []
    for n in nlist:
        print(n)
        i_init = np.argmin(abs(m.nrange/m.n0 - n))
        init = np.concatenate((m.rho[i_init, 2:], [m.mu_e[i_init]]))
        Es = []
        Ds = []
        Rs = []
        concs = []
        for f in frange:
            res, E = m.E_NF(n*m.n0, f, init)
            init = res[:-1].copy()
            n_n = n*m.n0 - sum(init[:-1])
            conc = np.insert(init[:-1], 0, n_n)
            concs.append(conc)
            Es.append(E)
#             Rs.append(init[-1]/n)
            D = eos.wrap_func_feq_rho(f, conc, init[-1], m.C)
#             print(D)
            Ds.append(D)
        Ess.append(Es)
        concss.append(concs)
        Dss.append(Ds)

    Ess = np.array(Ess)
    concss = np.array(concss)
    Dss = np.array(Dss)
    return [nlist, frange], [Ess, concss, Dss]

def get_Ef(m, n, frange=np.linspace(0, 1, 100)):
    if not hasattr(m, 'rho'):
        raise ValueError('Reset the model explicitly to provide initial values!')

    
    i_init = np.argmin(abs(m.nrange/m.n0 - n))
    init = np.concatenate((m.rho[i_init, 2:], [m.mu_e[i_init]]))
    Es = []
    concs = []
    for f in frange:
        res, E = m.E_NF(n*m.n0, f, init)
        init = res[:-1].copy()
        n_n = n*m.n0 - sum(init[:-1])
        conc = np.insert(init[:-1], 0, n_n)
        concs.append(conc)
        Es.append(E[0])
    
    Es = np.array(Es)
    concs = np.array(concs)
    return frange, Es, concs

def getBranchDict(dd, out):
    kind = 0
    if dd > 0:
        kind = 1
    return {
        'is_min' : kind,
        'n' : np.array([o[0] for o in out]),
        'f' : np.array([o[1] for o in out]),
        'conc' : np.array([o[2] for o in out]),
        'nc' : np.array([o[3] for o in out]),
        'E' : np.array([o[4] for o in out]),
        'D' : np.array([o[5] for o in out]),
        'dD' : np.array([o[6] for o in out]),
        'mu_e' : np.array([o[7] for o in out])
    }

def D_consistent(m, n, f, init, C, iterations=30):
    res, E = m.E_NF(n*m.n0, f, init, iterations=iterations)
    init = res[:-1].copy()
    n_n = n*m.n0 - sum(init[:-1])
    conc = np.insert(init[:-1], 0, n_n)
    D = eos.wrap_func_feq_rho(f, conc, init[-1], m.C)
    return D

def getBranches(m, paths, iterations=100):
    '''
        Returns dictionaries of branch properties for a given model and paths in (n, f) plane
    '''

    if not hasattr(m, 'rho'):
        raise ValueError('Reset the model explicitly to provide initial values!')
    branches = []
    for p in paths:
        out = []
        nlist, frange = p.vertices.transpose()
        flag = 0
        for n,f  in zip(tqdm(nlist), frange):
    #         print(n, f)
            if not flag:
                i_init = np.argmin(abs(m.nrange/m.n0 - n))
                init = np.concatenate((m.rho[i_init, 2:], [m.mu_e[i_init]]))
                flag = 1
            res, E = m.E_NF(n*m.n0, f, init, iterations=iterations)
            init = res[:-1].copy()
            mu_e = res[-2]
            n_n = n*m.n0 - sum(init[:-1])
            conc = np.insert(init[:-1], 0, n_n)
            D = eos.wrap_func_feq_rho(f, conc, init[-1], m.C)
            # dD = derivative(lambda z: eos.wrap_func_feq_rho(z, conc, init[-1], m.C), f, dx=1e-4)
            dD = derivative(lambda z: D_consistent(m, n, z, init, m.C, iterations=iterations),
                f, dx=1e-3, order=3)
            out.append([
                n,
                f,
                conc,
                res[-1],
                E[0],
                D,
                dD,
                mu_e
            ])

        ## Branch initial sign to compare with
        dd0 = out[0][6]
        # print(dd0)
        
        ## Add new branch when changing sign
        current = []
        for o in out:
            current += [o]
            if o[6] * dd0 > 0:
                # print(dd0, o[6])
                dd0 = o[6]
            else:
                # print('!!!!', dd0, o[6])
                branches.append(getBranchDict(dd0, current))
                dd0 = o[6]
                current = []
        branches.append(getBranchDict(dd0, current))
    return branches

def refineBranches(m, branches, iterations=300):
    new_branches = []
    for b in branches:
        out = []
        for i, n in enumerate(b['n']*m.n0):
            f = b['f'][i]
            init = b['conc'][i][1:] # shift for f and n_n
            init = np.append(init, b['mu_e'][i])
            # init = np.insert(init, 0, f)
            res, info = eos.stepE_rho_withf(n, init, np.array([f]), len(init)+2, 
                iterations, 0., m.C, 10+len(init))
            #inserting results
            init = res[:-1].copy()

            mu_e = res[-2]
            f = res[-3]
            
            n_n = n - sum(init[:-2])
            conc = np.insert(init[:-2], 0, n_n)
            n_in = np.insert(conc, 0, f)
            # print(len(n_in), n_in)
            # print(n, sum(n_in[1:]))
            nc = res[-1]
            E = m.Efull(rlist = [n_in], nclist=[nc], mu_list=[mu_e])
            D = 0
            dD = 0
            out.append([
                n,
                f,
                conc,
                res[-1],
                E[0],
                D,
                dD,
                mu_e]
            )
        new_branches.append(getBranchDict(-0.5 + b['is_min'], out))
    return new_branches


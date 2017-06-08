import eosWrap as eos
import numpy as np
import Models2
m = Models2.KVOR()

def Efull(n_in, rc2, C):
    mu_e = [eos.mu(r, 1, C) - eos.mu(r, 2, C) for r in n_in]
    E_free = np.array([eos.E_rho_free(n, r2, C) for n, r2 in zip(n_in, rc2)])
    n_e = [(mue**2 - m.m_e**2)**(1.5) / (3 * np.pi**2) if mue > m.m_e else 0 for mue in mu_e]
    n_mu = [(mue**2 - m.m_mu**2)**(1.5) / (3 * np.pi**2) if mue > m.m_mu else 0 for mue in mu_e]
    for i, ne in enumerate(n_e):
        E_free[i] += eos.kineticInt(ne, m.m_e, 2.)
        E_free[i] += eos.kineticInt(n_mu[i], m.m_mu, 2.)
    return np.array(E_free)

def solve(nrange, init, C, offset=None):
    rhos = []
    for _n in nrange:
        if offset:
            init += offset
        res = eos.step_free(_n, init, np.array([init[0]]), 3, 1500, C)
        print(_n, res)
        init = res.copy()
        rhos.append(res)

    rhos = np.array(rhos)

    n_in = []
    for i, r in enumerate(rhos):
        n_in.append(np.insert(r[:-1], 1, nrange[i] - r[1]))
    n_in = np.array(n_in)

    mu_e = [eos.mu(r, 1, C) - eos.mu(r, 2, C) for r in n_in]

    return n_in, rhos[:, -1], np.array(mu_e)

def get_Efr(frange, rcrange, init, C):
    n_0, rc2_0, f_0, n_p_0 = init
    E_rc2_f = []
    np_rc2_f = []

    for f in frange:
        n_p_init = n_p_0
        res = []
        res_np = []
        for rc2 in rcrange:
            n_p_init = eos.step_free_nfrc(n_0, n_p_init, rc2, f, 30, C)
            res_np.append(n_p_init)
            res.append(
                Efull([np.array([f, n_0 - n_p_init, n_p_init])], [rc2], C)[0])
        E_rc2_f.append(res)
        np_rc2_f.append(res_np)
    E_rc2_f = np.array(E_rc2_f)
    np_rc2_f = np.array(np_rc2_f)

    return E_rc2_f
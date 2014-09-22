import eosWrap as eos
import numpy as np
import multiprocessing as mp
import multi_progress as prog
from numpy import dtype
class Wrapper():
    def __init__(self, C):
        self.C = C
        self.m_pi = 135.0
        self.n0 = (197.33/self.m_pi)**3 * 0.16
        self.set = False
        self.driverSet = False
        self.const = 197.33*(0.16/self.n0)**(4.0/3.0)*self.m_pi**(-4.0)
    
    def reset(self, hyper=False, nmax = 4.0, npoints = 400, iter=30):
        if hyper:
            init = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        else:
            init = np.array([0.0], dtype=dtype('float'))
            
        self.n = np.linspace(0.01, nmax, npoints)
        last = init
        self.rho=[]

#         p = mp.Pool()
        writer = prog.Writer((0,0))
        pb = prog.ProgressBar(fd=writer)
        pb.start()
#         get_res = lambda z: self.rho = z
#         rs = p.map_async(lambda z: eos.stepE(z, last, len(last), iter, self.C),
#                 self.n, callback=get_res)
#         p.close()
#         while True:
#             if rs.ready():
#                 break
#             else:
#                 pb.update(100*rs._number_left/len(self.n))
#         
#         for i in rs:
#             self.rho.append(i)
        sp = 1 + self.C.sprime
        f = np.array([0.0 for i in range(sp)])
        for _i, i in enumerate(self.n):
            if abs(i - 3.53) < 1e-2:
                pass
            
            pb.update(100 * _i / len(self.n))
            rho = eos.stepE(i, last, f, len(last), iter, self.C)
            self.rho.append(rho.copy()) 
            self.rho[_i]=np.insert(self.rho[_i],0, i - np.sum(rho))
            f = eos.f_eq(self.rho[_i], f, sp, self.C)
            for j in range(sp):
                self.rho[_i]=np.insert(self.rho[_i], j, f[j])
            last = np.array(rho)
#             temp = self.rho[_i]
#             temp[3] = 0.0
#             print eos.mu(self.rho[_i], 1, self.C), eos.mu(temp, 3, self.C)
    
        self.rho=np.array(self.rho)
        _E = map(lambda z: eos.E(z, self.C), self.rho)
        dn = self.n[1] - self.n[0]
        self.P = self.n[1:]*np.diff(_E)/dn - _E[1:]
        self.E = np.ascontiguousarray(np.array(_E[1:]))
        self.n = np.ascontiguousarray(self.n[1:])
        self.rho = np.ascontiguousarray(self.rho[1:])
        self.set = True
    
    def EPN(self):
        if self.set:
            return self.E*self.m_pi**4, self.P*self.m_pi**4, self.n
        else:
            print 'Wrapper is not set!'
        
    def setDriver(self):
        if not self.set:
            print 'Wrapper is not set!'
            return 
        
        self.dr = eos.KVDriver()
        self.dr.set(self.E*self.m_pi**4, self.P*self.m_pi**4, self.n)
        self.driverSet = True
    
    def solve(self, f0=0.195, E0=-16.0, K0=275.0, J0=32.0, iter=1000, mu_scale = 1):
        sprime_old = self.C.sprime
        self.C.sprime = 0
        res = eos.solve(self.n0, E0, f0, K0, J0, self.C, iter, mu_scale)
        print res
        if mu_scale < 50000 and res == 5:
            self.solve(f0, E0, K0, J0, iter, 5*mu_scale)
        self.C.sprime = sprime_old 
    
    def stars(self, nmin = 0.5, nmax = 4.0, npoints=20):
        if not self.set:
            print 'Wrapper is not set!'
            return
        
        if not self.driverSet:
            print 'Driver is not set!'
            return 
        
        self.n_star = np.linspace(nmin, nmax, npoints)
        res = np.array(map(lambda z: eos.star(z, 2, self.dr), self.n_star))
        return self.n_star, res[:,0], res[:,1]
    
    def Esymm(self, n, f=0.0):
        res = []
        f = 0.0
        for z in n:
            f, = eos.f_eq(np.array([z/2,z/2]), np.array([f]), 1, self.C)
#             res.append((eos.E(np.array([f, z/2, z/2]), self.C))/z - self.C.M[0])
            res.append((eos.E(np.array([f, z/2, z/2]), self.C)))
        return np.array(res)
        
    def Psymm(self, n):
        E = self.Esymm(n)
        res = n[1:]*np.diff(E)/(n[1]-n[0]) - E[1:]
        return self.m_pi**4 * self.const*res
        
    
from scipy import interpolate, integrate
import numpy as np
from pylab import pi

class EnthalpySolver:
    
    def __init__(self, E=None, P=None, nhpoints=100):
        if E is not None and P is not None:
            self.setEos(E, P)
            self.nhpoints = nhpoints
        else:
            self.isSet = 0
        
        self.E_const = 4.898007e-4
        self.M_sun = 1.4766
        
    def equation(self, y, h):
        if not self.isSet:
            print "EoS is not set!"
            exit()
            
        r = y[0]
        m = y[1]
        if h < 0:
            h = 0.
#         print 'h=', h, 'r=',r,'m=',m
#         print 'P(h)=', self.PofH(h), 'E(h)=', self.EofH(h)
        dr = -(r * (r/self.M_sun - 2*m))/(m + r**3 * self.E_const * self.PofH(h))
        dm = - (self.E_const * self.EofH(h) *
                 r**3 * (r/self.M_sun - 2*m)) / (m + r**3 * 
                                              self.E_const * self.PofH(h))
        res = [dr, dm]

        return res 
    
       
    def setEos(self, E, P):
        self.E = E
        self.P = P
        self.EofP = interpolate.interp1d(P, E, kind='cubic')
        func_h = interpolate.interp1d(P, 1./(P + E), kind='cubic')
        self.h = lambda z: integrate.quad(func_h, P[0], z)[0]
        self.HofP = np.array(map(self.h, P))
        self.iHofP = interpolate.interp1d(P, self.HofP, kind='cubic', fill_value=self.HofP[-1])
        self.PofH = interpolate.interp1d(self.HofP, P, kind='cubic', fill_value=0.)
        self.EofH = interpolate.interp1d(self.HofP, E, kind='cubic', fill_value=0.)
        self.h_max = self.HofP[-1]
        self.isSet = 1
        
    def integrateOut(self, p_c):
        if not self.isSet:
            print "EoS is not set!"
            return -1

#         hc = self.HofP[-1]        
        hc = float(self.iHofP(p_c))

        print hc
#         exit()
        hrange = np.linspace(hc, 0., self.nhpoints)
        r0 = 1e-10
        MR = integrate.odeint(self.equation, [r0, 4*pi*r0**3*self.E_const*self.E[0]/3.],
                      hrange, rtol=1e-11, hmax=1e-3, tcrit=np.array([0., 0.]))
        return MR[:, 1], MR[:, 0], self.PofH(hrange), self.EofH(hrange)
    
    def integrateIn(self, m, r, hc):
        if not self.isSet:
            print "EoS is not set!"
            return -1

#         hc = self.HofP[-1]        
        print hc

        hrange = np.linspace(0., hc, 100*self.nhpoints)
        r0 = 1e-10
#         MR = integrate.odeint(self.equation, [r, m],
#                       hrange, rtol=1e-11, hmax=1e-3)
        ode = integrate.ode(lambda x, y: self.equation(y, x))
        ode.set_integrator('dopri5', atol=1e-8, rtol=1e-11, safety=0.2)
        ode.set_initial_value([r, m], 0.)
        dh = 1e-3
        MR = []
        while ode.successful() and ode.t < hc:
            ode.integrate(ode.t + dh)
            print ode.t, ode.y
            MR.append(ode.y)
        MR = np.array(MR)
        return MR[:, 1], MR[:, 0]
    
    def predictPE(self, m, r, e, p):
        D = 5.7207e-5
        e_new = e + 2.5 * (3 * m / (self.E_const * r**3) - e)
        p_new = (2.*pi/3) * D * (e + p) * (e + 3*p) * r**2
        p_new *= 1. + (2.*pi/3) * D * (4*e + 3*p) * r**2
        p_new += (pi/3) * D * (6*e + 11*p) * (3 * m / (self.E_const * r**3) - e) * r**2
        p_new += p
        return p_new, e_new
        
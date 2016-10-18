import numpy as np
from scipy.optimize import fsolve


class Sim():

    def __init__(self, q=1e-3, mdot=1e-6,alpha=.001,h=.05,nr=512,ri=1,ro=4,a=1,gamma=.5):



        self.q = q
        self.mdot = mdot
        self.alpha = alpha
        self.gamma = gamma
        self.h = h
        self.nr = nr
        self.ri = ri
        self.ro = ro
        self.a = a
        self.k = q**2/(self.alpha*self.h**5)
        kp = self.k*self.h**2
        self.gd = 1/(1+.04*self.k)
        self.gw = (.5*.5 + .16)*kp**(.25)


        self.r = np.exp(np.linspace(np.log(ri),np.log(ro),nr))

        self.dr = np.diff(np.log(self.r))[0]

        self.sig0 = self.mdot/(3*np.pi*self.nu(self.r))
        self.l = np.sqrt(self.r)

        self.xd = self.xdf(q)
        self.wd = self.wdf(q)

        self.dtr = self.dtrf(self.r)


    def nu(self,x):
        return self.alpha*self.h*self.h* x**(self.gamma)

    def xdf(self, q):
        k = q**2/(self.alpha*self.h**5)
        kp = k*self.h**2
        gd = 1/(1+.04*k)
        dr1 = (.25*gd + .08)*kp**(.25)
        dr2 = .33*kp**(.25)
        return (dr1 + dr2)/2
    def wdf(self,q):
        k = q**2/(self.alpha*self.h**5)
        kp = k*self.h**2
        gd = 1/(1+.04*k)
        dr1 = (.25*gd + .08)*kp**(.25)
        dr2 = .33*kp**(.25)
        return (dr2-dr1)

    def dtrf(self,x,fnl=.4):
        delta=  self.h*self.a
        res = (.4*fnl * 2*np.pi * (1./(x/self.a-1))**4)*np.array(abs(x-self.a)>=1.3*delta).astype(int)
        res[np.isnan(res)] = 0
        return res

    def dep_func(self,x,q):
        dist = np.array(x/self.a-1)
        return ((dist>=self.xd-self.wd/2)&(dist<=self.xd+self.wd/2)).astype(int)* 1./self.wd
    def dep_func_int(self,x,q):
        dist = np.array(x/self.a-1)
        return ((dist>=self.xd-self.wd/2)&(dist<=self.xd+self.wd/2)).astype(int)* (dist -(self.xd-self.wd/2))/self.wd + (dist>self.xd+self.wd/2).astype(int)

    def func(self,T):
        dens = self.sig0*(1 + T*self.dep_func_int(self.r,self.q)/(self.mdot*self.l) - (self.l[0]/self.l)*(.04*self.k*self.gd))
        integ = (self.dr*self.dtr*dens).sum()
        q = np.sqrt(T/integ)

        return self.q/q - 1

    def solve(self):

        if self.q < self.h**3:
            fnl = 1
        else:
            fnl = .4

        T0 = 2*np.pi*fnl*self.a*self.q**2 / self.h**3
        T0 *= self.mdot/(3*np.pi*self.nu(self.a))

        T = fsolve(self.func,T0)

        self.dens_ans =  self.sig0*(1 + T*self.dep_func_int(self.r,self.q)/(self.mdot*self.l)- (self.l[0]/self.l)*(.04*self.k*self.gd))


        return self.dens_ans,T







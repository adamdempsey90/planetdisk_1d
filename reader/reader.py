import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz
import h5py
from collections import OrderedDict
from subprocess import call

class Parameters():
    key_vals=OrderedDict([
        ('nr',int),
        ('ri',float),
        ('ro',float),
        ('alpha',float),
        ('gamma',float),
        ('h',float),
        ('bc_lam_inner',float),
        ('bc_lam_outer',float),
        ('bc_mdot',float),
        ('flux_bc',bool),
        ('dt',float),
        ('cfl',float),
        ('nvisc',float),
        ('nt',int),
        ('release_time',float),
        ('start_ss',bool),
        ('read_initial_conditions',bool),
        ('planet_torque',bool),
        ('explicit_stepper',bool),
        ('move_planet',bool),
        ('move_planet_implicit',bool),
        ('gaussian',bool),
        ('symmetric_torque',bool),
        ('nonlocal_torque',bool),
        ('shock_dep',bool),
        ('forced_torque',bool),
        ('hs_visc',bool),
        ('one_sided',float),
        ('a',float),
        ('mp',float),
        ('G1',float),
        ('beta',float),
        ('delta',float),
        ('c',float),
        ('eps', float),
        ('xd',float),
        ('outputname',str)])

    comment_lines = {'nr':'#Grid Parameters:', 'alpha':'\n\n#Disk Parameters:','bc_lam_inner':'\n\n#Boundary Conditions:', 'dt':'\n\n#Time Parameters:','planet_torque':'\n\n#Planet Properties:'}

    def __init__(self,dataspace,fname='results.hdf5'):
           # We read in the from an hdf5 file
        setattr(self,'outputname',fname)
        for key,datatype in list(self.key_vals.items())[:-1]:
            try:
                setattr(self,key,dataspace[key].astype(datatype)[0])
            except AttributeError:
                setattr(self,key,getattr(dataspace,key))

    def dump_params(self,fname='params_py.in',**kargs):
        outname = kargs.pop('outputname',self.outputname)
        lines =[]
        for key in self.key_vals.keys():
            try:
                lines.append(self.comment_lines[key])
            except KeyError:
                pass

            if key in kargs.keys():
                val = kargs[key]
            else:
                val = getattr(self,key)
                if key == 'a':
                    val = self.at[0]
                if key == 'outputname':
                    val = outname

            lines.append(str(key) + ' = ' + str(val))

        lines = '\n'.join(lines) + '\n'

        with open(fname,'w') as f:
            f.write(lines)

    def run(self,fname='params_py.in',**kargs):
        outname = kargs.pop('outputname',self.outputname)

        self.dump_params(fname,outputname=outname,**kargs)
        call(['../pdisk1d',fname])
        return outname

class Sim(Parameters):
    def __init__(self,fname = 'results.hdf5'):

        f = h5py.File(fname,'r')

        Parameters.__init__(self,f['/Migration/Parameters']['Parameters'],fname=fname)

        sols = f['/Migration/Solution']
        mesh = f['/Migration/Mesh']
        mat = f['/Migration/Matrix']
        ss = f['/Migration/SteadyState']
        self.t = sols['times'][:]
        self.at =sols['avals'][:]
        self.vs = sols['vs'][:]
        self.lam = sols['lam'][:].transpose()
        self.mdot = sols['mdot'][:].transpose()
        self.torque= sols['torque'][:].transpose()
        self.dTr = ss['dTr'][:].transpose()

        self.rc = mesh['rc'][:]
        self.l = np.sqrt(self.rc)
        self.dr = mesh['dr'][:]
        self.tauc = mesh['tauc'][:]
        self.taum = mesh['taumin'][:]
        self.dlr = self.dr/self.rc
        self.rmin=mesh['rmin'][:]
        self.lam0 = mesh['lami'][:]
        self.mdot0 = mesh['mdoti'][:]
        self.nu_grid = mesh['nu_grid'][:]
        self.dep_func = mesh['dep_func'][:]
        self.main_diag = mat['md'][:]
        self.lower_diag = mat['ld'][:]
        self.upper_diag = mat['ud'][:]
        self.rhs = mat['fm'][:]
        self.u = mat['u'][:]
        self.w = mat['w'][:]

        self.lamex = self.dTr * self.lam
        self.ilamex = self.lamex * self.dr[:,np.newaxis]
        self.ilamdep = self.torque * self.dr[:,np.newaxis]

        try:
            self.col = mat['col'][:]
        except KeyError:
            print('Could not find col entry in file')

        self.lamnp = ss['lam0'][:].transpose()
        self.lam_ss = ss['lam_ss'][:].transpose()
        self.lamp = ss['lamp'][:].transpose()
        self.mdot_ss = ss['mdot_ss'][:]
        self.ivals_ss = ss['ivals_ss'][:].transpose()
        self.kvals_ss = ss['kvals_ss'][:].transpose()

        self.vs_ss = ss['vs_ss'][:]
        self.eff = ss['eff'][:]
 #       self.mdot0 = self.mdot_ss/self.eff
        f.close()
        self.fac = 2*np.pi*self.rc
        self.sigma = self.lam/self.fac[:,np.newaxis]
        self.fh = self.lam*(1.5*self.nu(self.rc)/np.sqrt(self.rc) )[:,np.newaxis]
#
#        self.t = dat_p[:,0]
#        self.nt = len(self.t)
#        self.at = dat_p[:,1]
#        self.vst = dat_p[:,2]
#
#        dat = np.loadtxt(fname)
#        self.nr = dat[0,0]
#        self.rc = dat[0,1:]
#        self.rm = dat[1,1:]
        self.mth = self.h**3
        self.q = self.mp*self.mth
        self.rh = (self.q/3.)**(1./3) * self.a
        self.disk_mass = np.dot(self.dr,self.lam0)
#        self.dr = dat[2,1:]
        self.tvisc = self.ro**2/(self.nu(self.ro))
        self.tviscp = self.at**2/(self.nu(self.at))
        self.lami = self.bc_lam_inner
        self.lamo = self.bc_lam_outer
#        self.K = self.mp * self.h/self.alpha
#        self.B = 2*self.at*(self.lamo-self.lami)/(self.mp*self.mth)
#        self.Bfac = (np.sqrt(self.ro)-np.sqrt(self.ri))/(np.sqrt(self.at)-np.sqrt(self.ri))
        self.Ap = self.lamp[-1,:]*np.sqrt(self.ro)
        self.vp_visc = -1.5*self.nu(self.at)/self.at
        self.A = (self.torque*self.lam)
        for i in range(self.lam.shape[1]):
            self.A[:,i] *= self.dlr/self.bc_mdot
        self.A = np.trapz(self.A,axis=0)

        self.K = self.eps*(self.mp*self.mth)**2/(self.alpha*self.h**5)
        self.B = 2 * self.at * self.bc_mdot/(-self.vr_nu(self.at))
        self.B /= self.q
        self.freduc = self.B * self.Ap/np.sqrt(self.at)

        self.Fh = 1.5*self.lam*( (self.nu(self.rc)/np.sqrt(self.rc))[:,np.newaxis])
        self.dFh = self.grad(self.Fh)*2*np.sqrt(self.rc[:,np.newaxis])
        self.Fh_ss = 1.5*self.lam_ss*( (self.nu(self.rc)/np.sqrt(self.rc))[:,np.newaxis])
        self.dFh_ss = self.grad(self.Fh_ss)*2*np.sqrt(self.rc[:,np.newaxis])
        self.mdotl = np.zeros(self.mdot.shape)
        self.mdotl[0,:] = (self.mdot[0,:]*self.l[0])
        self.mdotl[1:,:] = cumtrapz(self.mdot*(np.sqrt(self.rc)*.5*self.dlr)[:,np.newaxis],axis=0)
        for i in range(1,self.mdot.shape[0]):
            self.mdotl[i,:] += self.mdotl[0,:]
        self.iT = np.zeros(self.mdot.shape)
        self.iT[1:,:] = cumtrapz(self.torque*self.lam*self.dlr[:,np.newaxis],axis=0)
#        self.dTr = dat[3,1:]
#        self.lami = dat[4,0]
#        self.lam0 = dat[4,1:]
#        self.lamo = dat[5,0]
#        self.mdot0 = dat[5,1:]
#        self.lams = dat[6:-self.nt,1:].transpose()
#        self.mdots = dat[-self.nt:,1:].transpose()
        self.vr = -self.mdot/self.lam
        self.vr0 = np.zeros(self.vr.shape)
        self.vr_ss = np.zeros(self.vr.shape)
#        self.lamp = np.zeros(self.lam.shape)
        for i in range(self.vr0.shape[1]):
#            self.lamp[:,i] = self.lam[:,i]/self.lam0
            self.vr0[:,i] = -self.mdot0[i]/self.lam0
            self.vr_ss[:,i] = -self.mdot_ss[i]/self.lam_ss[:,i]
        self.beta_reduc = self.vr_ss/self.vr0

        self.set_ss()
#        self.vr0 = -self.mdot0/self.lam0
#        self.lamp = np.zeros(self.lam.shape)
#        for i in range(self.lam.shape[1]):
#            self.lamp[:,i] = (self.lam[:,i]-self.lam0)/self.lam0
    def set_ss(self):
        if self.flux_bc:
            mdot = self.bc_mdot

            for i in range(self.lam_ss.shape[1]):
                self.vr0[:,i] = self.vr_nu(self.rc) #/(1-np.sqrt(self.ri/self.rc))
                #self.lam_ss[:,i] /= self.eff[i]

         #   self.lam0 = -mdot/self.vr_nu(self.rc)
        else:
            for i in range(self.lam_ss.shape[1]):
                self.vr0[:,i] = self.vr_nu(self.rc)/(1-np.sqrt(self.ri/self.rc))
                self.lam_ss[:,i]  /= self.eff[i]
         #   self.lam0 = -mdot/self.vr_nu(self.rc)
    def nu(self,x):
        return self.alpha*self.h*self.h * pow(x,self.gamma)
    def vr_nu(self,x):
#        return -1.5 * self.nu(x)/x /(1 - np.sqrt(self.ri/x))
        return -1.5*self.nu(x)/x
    def grad(self,q,axis=0):
        res = np.zeros(q.shape)

        if axis == 0:
            res[1:-1,:] = (q[2:,:] - q[:-2,:])/(self.rc[2:]-self.rc[:-2])[:,np.newaxis]
            res[0,:] = (q[1,:]-q[0,:])/(self.rc[1]-self.rc[0])
            res[-1,:] = (q[-1,:]-q[-2,:])/(self.rc[-1]-self.rc[-2])
        elif axis == 1:
            res[:,1:-1] = (q[:,2:] - q[:,:-2])/(self.rc[2:]-self.rc[:-2])[:,np.newaxis]
            res[:,0] = (q[:,1]-q[:,0])/(self.rc[1]-self.rc[0])
            res[:,-1] = (q[:,-1]-q[:,-2])/(self.rc[-1]-self.rc[-2])
        return res

    def plot(self,q='lam',i=-1,logx=False,logy=False,norm=None,ax=None,fig=None):

        try:
            dat = getattr(self,q)[:,i].copy()
        except AttributeError:
            print('Q not found!')
            return
        except IndexError:
            print('{:d} not a valid index!'.format(i))
            return

        if ax is None:
            fig= plt.figure()
            ax = fig.add_subplot(111)

        if norm is not None:
            dat /= norm

        ax.plot(self.rc,dat,'-k')

        if logx:
            ax.set_xscale('log')
        if logy:
            ax.set_yscale('log')

        ax.set_xlabel('Radius [AU]',fontsize=15)




    def animate(self,tend,skip,tstart=0,q='lam',logx = True,logy=True,ylims=None):
        fig=plt.figure()
        ax=fig.add_subplot(111)
        inds = (self.t <= tend)&(self.t >= tstart)
        ax.set_xlabel('$r$',fontsize=20);

        if q == 'lam':
            ax.set_ylabel('$\\lambda$',fontsize=20)
            line, = ax.plot(self.rc,self.lam0)
            linep, = ax.plot(self.at[0],self.lam0[self.rc>=self.at[0]][0],'o',markersize=10)
            ax.plot(self.rc,self.lam_ss[:,-1],'--r')
            ax.plot(self.rc,self.lam0,'--k')
            dat = self.lam[:,inds]
            dat = dat[:,::skip]
            if logy:
                ax.set_yscale('log')
        elif q == 'mdot':
            ax.set_ylabel('$\\dot{M}$',fontsize=20)
            line, = ax.plot(self.rc,self.mdot[:,0]/self.bc_mdot)
            linep, = ax.plot(self.at[0],1,'o',markersize=10)

            #ax.axhline(1,color='k',linestyle='--')
            dat = self.mdot[:,inds]/self.bc_mdot
            dat = dat[:,::skip]
        elif q == 'torque':
            ax.set_ylabel('$ \\Lambda(r)$',fontsize=20)
            line, = ax.plot(self.rc,self.torque[:,0])
            linep, = ax.plot(self.at[0],self.torque[:,0][self.rc>=self.at[0]][0],'o',markersize=10)
            dat = self.torque[:,inds][:,::skip]

        elif q == 'sigma':
            fac = 2*np.pi*self.rc
            ax.set_ylabel('$ \\Sigma(r)$',fontsize=20)
            line, = ax.plot(self.rc,self.lam0/fac)
            linep, = ax.plot(self.at[0],(self.lam0/fac)[self.rc>=self.at[0]][0],'o',markersize=10)
            liness,=ax.plot(self.rc,self.lam_ss[:,0]/fac,'--r')
            ax.plot(self.rc,self.lam0/fac,'--k')
            ax.plot(self.rc,self.lamnp[:,0]/fac,'-k')
            dat = self.lam[:,inds]
            dat = dat[:,::skip]
            for i in range(dat.shape[1]):
                dat[:,i] /= fac

            if logy:
                ax.set_yscale('log')
        elif q == 'lamp':
            ax.set_ylabel('$ \\Sigma/\\Sigma_0$',fontsize=20)
            line,=  ax.plot(self.rc,self.lamp[:,0]+1)
            linep, = ax.plot(self.at[0],(self.lamp[:,0]+1)[self.rc>=self.at[0]][0],'o',markersize=10)
            liness, = ax.plot(self.rc,self.lam_ss[:,0]/self.lamnp[:,0],'--r')
            ax.plot(self.rc,self.lamp[:,0]+1,'--k')
            dat = self.lamp[:,inds][:,::skip] + 1

        elif q == 'mass':
            ax.set_ylabel('Ring Mass',fontsize=15)
            line, = ax.plot(self.rc,self.lam0*self.dr)
            linep, = ax.plot(self.at[0],(self.lam0*self.dr)[self.rc>=self.at[0]][0],'o',markersize=10)
            liness,=ax.plot(self.rc,self.lam_ss[:,0]*self.dr,'--r')
            ax.plot(self.rc,self.lam_ss[:,-1]*self.dr,'--r')
            ax.plot(self.rc,self.lam0*self.dr,'--k')
            dat = self.lam[:,inds]*self.dr[:,np.newaxis]
            dat = dat[:,::skip]
            if logy:
                ax.set_yscale('log')
        elif q == 'cmass':
            ax.set_ylabel('Cumulative Interior Mass',fontsize=15)
            line, = ax.plot(self.rc,np.cumsum(self.lam0*self.dr))
            linep, = ax.plot(self.at[0],np.cumsum(self.lam0*self.dr)[self.rc>=self.at[0]][0],'o',markersize=10)
            liness,=ax.plot(self.rc,np.cumsum(self.lam_ss[:,0]*self.dr),'--r')
            ax.plot(self.rc,np.cumsum(self.lam_ss[:,-1]*self.dr),'--r')
            ax.plot(self.rc,np.cumsum(self.lam0*self.dr),'--k')
            dat = np.cumsum(self.lam[:,inds]*self.dr[:,np.newaxis],axis=0)
            dat = dat[:,::skip]
            if logy:
                ax.set_yscale('log')
        else:
            print('q={} is not a valid option'.format(q))
            return

        if logx:
            ax.set_xscale('log')
        ax.set_ylim((dat.min(),dat.max()))

        times = self.t[inds][::skip]

        avals = self.at[inds][::skip]
        if ylims != None:
            ax.set_ylim(ylims)
        for i,t in enumerate(times):
            line.set_ydata(dat[:,i])
            linep.set_xdata(avals[i])
            linep.set_ydata(dat[:,i][self.rc>=avals[i]][0])
            if q == 'sigma':
                liness.set_ydata(self.lam_ss[:,inds][:,::skip][:,i]/fac)
            if q == 'mass':
                liness.set_ydata(self.lam_ss[:,inds][:,::skip][:,i]*self.dr)
            if q == 'cmass':
                liness.set_ydata(np.cumsum(self.lam_ss[:,inds][:,::skip][:,i]*self.dr))
            ax.set_title('t = {:.2e} viscous times'.format(t/self.tvisc))
            plt.pause(1e-4)
            plt.draw()


    def time_series(self,axes=None,fig=None,scale=True):
        if fig == None:
            fig,axes = plt.subplots(2,3,sharex='col')
        fig.subplots_adjust(hspace=0)
        axes[0,0].semilogx(self.t,self.at,'.-')

        if scale:
            dat = self.vs / self.vp_visc
            vlabel = '$\\dot{a}/v_r^{visc}$'
        else:
            dat = self.vs
            vlabel = '$\\dot{a}$'

        axes[1,0].semilogx(self.t,dat,'.-')

        axes[1,0].set_yscale('log')
        axes[0,0].set_ylabel('$a$',fontsize=20)
        axes[1,0].set_ylabel(vlabel,fontsize=20)
        axes[1,0].set_xlabel('$t$',fontsize=20)

        if (self.at >= self.rc[-1]).any():
            axes[0,0].axhline(self.rc[-1],color='k')
        if (self.at <= self.rc[0]).any():
            axes[0,0].axhline(self.rc[0],color='k')

        axes[1,1].plot(self.at,dat,'.-')
        axes[1,2].plot(self.B,dat,'.-')
        axes[1,1].set_xlabel('$a$',fontsize=20)
        axes[1,1].set_ylabel(vlabel,fontsize=20)
        if scale:
            axes[1,1].set_yscale('log')
            axes[1,2].set_yscale('log')
            axes[1,1].plot(self.at,self.freduc,'+r',label='steady state pred.') #self.freduc,'+r')
            axes[1,2].plot(self.B,self.freduc,'+r')
        else:
            axes[1,1].plot(self.at,self.vs_ss,'+r',label='steady state pred.') #self.freduc*self.vr_nu(self.at),'+r')
            axes[1,2].plot(self.B,self.vs_ss,'+r')

        axes[1,1].legend(loc='best')
        axes[0,1].plot(self.at,self.A,'.-',label='A')
        axes[0,1].plot(self.at,self.B,'.-',label='B')
        axes[0,1].legend(loc='best')
        axes[0,1].set_yscale('log')


        axes[0,2].plot(self.B,self.A/np.sqrt(self.at),'.-')
        axes[0,2].plot(self.B,self.Ap/np.sqrt(self.at),'+r')
        axes[1,2].set_xlabel('$B$',fontsize=20)
        axes[0,2].set_ylabel('$A/\\sqrt{a}$',fontsize=20)


        add_line = lambda ax,val,orient,c,w: ax.axhline(val,color=c,linewidth=w) if orient == 'y' else ax.axvline(val,color=c,linewidth=w)


        lowx, highx = axes[0,0].get_xlim()
        lowy, highy = axes[0,0].get_ylim()
        if lowy <= self.ri <= highy:
            add_line(axes[0,0],self.ri,'y','k',2)

        lowx, highx = axes[0,1].get_xlim()
        lowy, highy = axes[0,1].get_ylim()
        if lowy <= 1 <= highy:
            add_line(axes[0,1],1,'y','k',2)
        if lowx <= self.ri <= highx:
            add_line(axes[0,1],self.ri,'x','k',2)

        lowx, highx = axes[1,0].get_xlim()
        lowy, highy = axes[1,0].get_ylim()
        if lowy <= 1 <= highy:
            add_line(axes[1,0],1,'y','k',2)

        lowx, highx = axes[1,1].get_xlim()
        lowy, highy = axes[1,1].get_ylim()
        if lowy <= 1 <= highy:
            add_line(axes[1,1],1,'y','k',2)
        if lowx <= self.ri <= highx:
            add_line(axes[1,1],self.ri,'x','k',2)

        lowx, highx = axes[1,2].get_xlim()
        lowy, highy = axes[1,2].get_ylim()
        if lowy <= 1 <= highy:
            add_line(axes[1,2],1,'y','k',2)
        if lowx <= 1 <= highx:
            add_line(axes[1,2],1,'x','k',2)

        lowx, highx = axes[0,2].get_xlim()
        lowy, highy = axes[0,2].get_ylim()
        if lowx <= 1 <= highx:
            add_line(axes[0,2],1,'x','k',2)

        for ax in axes.flatten():
            ax.minorticks_on()
            ax.grid()

        fig.canvas.draw()
        return axes,fig
    def write_lam_to_file(self,r,lam,interpolate=False):
        if interpolate:
            lam = interp1d(r,lam)(self.rc)
            r = self.rc

        with open("lambda_init.dat","w") as f:
            lines = ["{:.12f}\t{:.12f}".format(x,l) for x,l in zip(r,lam)]
            f.write('\n'.join(lines))

    def load_steadystate(self,directory=''):
        if len(directory)>0 and directory[-1] != '/':
            directory += '/'
        dat = np.loadtxt(directory+'results.dat')
        self.ss_r = dat[:,0]
        self.ss_lam = dat[:,1]
        self.ss_vr = dat[:,2]
        self.ss_dTr = dat[:,3]
        self.ss_mdot = -dat[10,1]*dat[10,2]

def plaw(x,a,b):
    return a*x**b


class Field():
    def __init__(self,fname='outputs/results.hdf5'):
        dat = h5py.File(fname)
        sols = dat['/Migration/Solution']
        mesh = dat['/Migration/Mesh']
        mat = dat['/Migration/Matrix']
        self.rc = mesh['rc'][:]
        self.rmin = mesh['rmin'][:]
        self.nr = len(self.rc)
        self.dr = mesh['dr'][:]
        self.dep_func = mesh['dep_func'][:].transpose()
        self.t =  sols['times'][:]
        self.mdot = sols['mdot'][:].transpose()
        self.torque = sols['torque'][:].T
        self.mdL = sols['mdL'][:].T
        self.mdR = sols['mdR'][:].T
        self.lam = sols['lam'][:].T
        self.dTr = sols['dTr'][:].T * self.lam
        self.sigma = self.lam/(2*np.pi*self.rc[:,np.newaxis])
        self.drlam = np.zeros(self.lam.shape)
        self.drlam[1:-1,:] = (self.lam[2:,:]-self.lam[:-2,:])/(self.rc[2:]-self.rc[:-2])[:,np.newaxis]
        self.drlam[0,:] = (self.lam[1,:]-self.lam[0,:])/(self.rc[1]-self.rc[0])
        self.drlam[-1,:] = (self.lam[-1,:]-self.lam[-2,:])/(self.rc[-1]-self.rc[-2])
        self.dtlam = np.diff(self.lam,axis=1)/np.diff(self.t)
        self.mdot_c = np.zeros((self.nr,len(self.t)-1))
        for i in range(len(self.t)-1):
            self.mdot_c[:,i] = ((self.dtlam[:,i]*self.dr)[::-1].cumsum()+self.mdot[-1,i])[::-1]
        self.md = mat['md'][:]
        self.ld = mat['ld'][:]
        self.ud = mat['ud'][:]
        self.fm = mat['fm'][:]


        temp = mat['u'][:]
        self.uL = temp[:self.nr]
        self.uR = temp[-self.nr:]
        temp = mat['w'][:]
        self.wL = temp[:self.nr]
        self.wR = temp[-self.nr:]
        dat.close()
    def plotall(self,q='sigma',ax=None,norm=None,logx=False,logy=False,**kargs):
        try:
            dat = getattr(self,q).copy()
        except AttributeError:
            print('{} is ot a valid field!'.format(q))
            return
        if norm is not None:
            try:
                dat /= norm[:,np.newaxis]
            except (TypeError, ValueError) as e:
                dat /= norm

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        ax.plot(self.rc,dat,**kargs)
        ax.set_xlabel('Radius [AU]',fontsize=15)
        if logx:
            ax.set_xscale('log')
        if logy:
            ax.set_yscale('log')

    def loglog(self,**kargs):
        self.plot(logx=True,logy=True,**kargs)
    def semilogx(self,**kargs):
        self.plot(logx=True,logy=False,**kargs)
    def semilogy(self,**kargs):
        self.plot(logx=False,logy=True,**kargs)
    def plot(self,q='sigma',i=-1,ax=None,initial=False,norm=None,logy=False,logx=False,**kargs):
        try:
            dat = getattr(self,q)[:,i].copy()
            if initial:
                dat0 = getattr(self,q)[:,0].copy()
        except AttributeError:
            print('{} is ot a valid field!'.format(q))
            return
        if norm is not None:
            try:
                dat /= norm[:,np.newaxis]
                if initial:
                    dat0 /= norm[:,np.newaxis]
            except (TypeError,ValueError) as e:
                dat /= norm
                if initial:
                    dat0 /= norm
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        ax.plot(self.rc,dat,**kargs)
        if initial:
            ax.plot(self.rc,dat0,'--k')
        ax.set_xlabel('Radius [AU]',fontsize=15)
        if logx:
            ax.set_xscale('log')
        if logy:
            ax.set_yscale('log')
    def plot_itorque(self,ax=None,logx=False,i=-1,a=1):
        if ax is None:
            fig = plt.figure()
            ax=fig.add_subplot(111)


        TeR = (self.dr*self.dTr[:,i])[self.rc>=a].cumsum()
        TeL = -(self.dr*self.dTr[:,i])[self.rc<=a][::-1].cumsum()
        TdR = (self.dr*self.torque[:,i])[self.rc>=a].cumsum()
        TdL = -(self.dr*self.torque[:,i])[self.rc<=a][::-1].cumsum()

        ax.plot(self.rc[self.rc>=a],TeR,'-k')
        ax.plot(self.rc[self.rc<=a][::-1],TeL,'-k')
        ax.plot(self.rc[self.rc>=a],TdR,'-r')
        ax.plot(self.rc[self.rc<=a][::-1],TdL,'-r')

        ax.set_xlabel('Radius [AU]',fontsize=15)
    def zero_torque(self,ax=None,remove=False,fac=3*np.pi*.001*.05*.05,gamma=0):
        md = self.mdot_c[:,-1]*np.sqrt(self.rc) - (fac*self.rc**gamma*self.sigma[:,-1])
        if remove:
            md = self.mdot_c[0,-1]*np.sqrt(self.rc[0]) - (fac*self.rc[0]**gamma *self.sigma[0,-1])
        dt = -(self.dr*self.torque[:,-1]).cumsum()
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        ax.plot(self.rc,md,'-b')
        ax.plot(self.rc,dt,'-r')
        ax.set_xlabel('Radius [AU]',fontsize=15)


    def summary(self,axes=None,logx=False,logy=True):
        if axes is None:
            fig,axes= plt.subplots(2,1,sharex=True)

        axes[0].plot(self.rc,self.sigma[:,-1]/self.sigma[:,0])
        axes[1].plot(self.rc,self.mdot_c[:,-1])

        for ax in axes:
            if logx:
                ax.set_xscale('log')
            ax.minorticks_on()

        if logy:
            axes[0].set_yscale('log')

        axes[1].set_xlabel('Radius [AU]',fontsize=15)

        axes[1].set_ylabel('$\\dot{M}$',fontsize=20)
        axes[0].set_ylabel('$\\Sigma/\\Sigma_0$',fontsize=20)
    def evaulate_matrix(self):
        fac1 = np.zeros(self.rc.shape)

        for i in range(1,self.nr-1):
            fac1[i] = self.md[i]*self.lam[i,-1] + self.ld[i-1]*self.lam[i-1,-1]+self.ud[i]*self.lam[i+1,-1]

        fac1[0] = self.lam[0,-1]*self.md[0] + self.lam[1,-1]*self.ud[0]
        fac1[-1] = self.lam[-1,-1]*self.md[-1] + self.lam[-2,-1]*self.ld[-1]

        fac2 = self.uL *np.dot(self.wL,self.lam[:,-1]) + self.uR*np.dot(self.wR,self.lam[:,-1])

        plt.figure()
        plt.plot(self.rc,fac1,label='main')
        plt.plot(self.rc,fac2,label='uw')
        plt.plot(self.rc,self.fm,label='F')
        plt.legend()
        plt.figure()
        plt.plot(self.rc,fac1+fac2-self.fm ,'-b')

def gaussian_ring(r,t, alpha=.01,h=.05,gamma=.5,r1=1.,m=1.):
    from scipy.special import iv
    ts = 12*alpha*h*h*r1**gamma * t/r1**2
    norm = m/(np.pi*r1**2)
    u = r/r1
    return norm * iv(.25,2*u/ts) * np.exp(-(1 + u*u)/ts) / (u**(.25) * ts)

def self_similar(sigma0, r,t, alpha=.01,h=.05,gamma=.5,r1=1.):
    ts = 1./(3*(2-gamma)**2) * r1**2/(alpha*h**2*r1**gamma)
    tval = t/ts + 1
    return sigma0 * tval**( - (2.5-gamma)/(2-gamma)) * np.exp(- (1-tval)/tval * r/r1)

def self_sim_bc(gamma, rin,rout):
    res = np.array([rin,rout])
    return (gamma-1) + (2-gamma)*res**(2-gamma)



class Results():
    def __init__(self,fname='output_test.dat'):
        with open(fname,"r") as f:
            n,m,mstart,mend=np.fromfile(f,dtype=float,count=4).astype(int)
            self.TL = np.fromfile(f,dtype=float,count=m)
            self.TR = np.fromfile(f,dtype=float,count=m)
            self.r = np.fromfile(f,dtype=float,count=n)
            self.lamex = np.fromfile(f,dtype=float,count=n*m).reshape(m,n).T
            self.lamdep = np.fromfile(f,dtype=float,count=n*m).reshape(m,n).T
            self.drfw = np.fromfile(f,dtype=float,count=n*m).reshape(m,n).T
            self.fw = np.fromfile(f,dtype=float,count=n*m).reshape(m,n).T
            self.dppot = np.fromfile(f,dtype=float,count=n*m).reshape(m,n).T
            self.u= np.fromfile(f,dtype=complex,count=n*m).reshape(m,n).T
            self.v = np.fromfile(f,dtype=complex,count=n*m).reshape(m,n).T
            self.sig = np.fromfile(f,dtype=complex,count=n*m).reshape(m,n).T
            self.dbar = np.fromfile(f,dtype=float,count=n)
            self.n,self.nm,self.mstart,self.mend = n,m,mstart,mend
    def mode_plot(self,m):
        fig,axes=plt.subplots(1,3,figsize=(15,5))
        axes[0].plot(self.r,self.u.real[:,m-1],self.r,self.u.imag[:,m-1])
        axes[1].plot(self.r,self.v.real[:,m-1],self.r,self.v.imag[:,m-1])
        axes[2].plot(self.r,self.sig.real[:,m-1],self.r,self.sig.imag[:,m-1])
    def torque(self,m,tot=False):
        fig=plt.figure()
        ax = fig.add_subplot(111)
        if tot:
            ax.plot(self.r,self.lamex.sum(axis=1),'-k')
            ax.plot(self.r,self.drfw.sum(axis=1),'-r')
            ax.plot(self.r,self.lamdep.sum(axis=1),'-b')
        else:
            ax.plot(self.r,self.lamex[:,m-1],'-k')
            ax.plot(self.r,self.drfw[:,m-1],'-r')
            ax.plot(self.r,self.lamdep[:,m-1],'-b')


#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

def load_data():
    a = np.loadtxt('log.out')
    t_sho = a[:,0]  # Trap time
    tv_scalar = a[:,4]  # Scalar time for some other sims
    #tv_scalar = a[:,3] # Scalar time

    tv = tv_scalar
    d = np.loadtxt('trap.dat')
    xv = d[:,0]
    w = d[:,2]
    
    nt = tv.shape[0]
    nx = xv.size
    
    f = np.fromfile('fields.bin').reshape((nt,2,2,nx))

    return xv, t_sho, w, f

def plot_rho_tot_deformation(xv,tv,psi):
    rho = np.sum(psi**2,axis=-2)
    rhot = 0.5*(rho[:,1]+rho[:,0])

    rho_rel = (rhot-rhot[0,:]) # Fix this to normalize to initial density
    rho_dist = max(np.abs(np.min(rho_rel)), np.max(rho_rel))
    levels = np.linspace(-rho_dist,rho_dist,20,endpoint=True)
    
    f,a = plt.subplots()
    cnt = a.contourf(xv, tv, rho_rel, levels, cmap='RdBu')
    for c_ in cnt.collections:
        c_.set_edgecolor("face")
    cv = f.colorbar(cnt)

    a.set_xlim(-5.,5.)
    return f,a

def plot_rho_tot_max(w,tv,psi):
    rho = np.sum(psi**2,axis=-2)
    rhot = 0.5*(rho[:,1]+rho[:,0])
    rhot = rhot - rhot[0]
    
    rhot_s = np.sum(w[np.newaxis,:]*np.abs(rhot),axis=-1)
    
    f,a = plt.subplots()
    a.plot(tv,rhot_s)
    a.plot(tv,np.max(np.abs(rhot),axis=-1))
    return f,a

def plot_rho_rel(xv,tv,psi,nu):
    f,a = plt.subplots()

    rho = np.sum(psi**2,axis=-2)
    drho = 0.5*(rho[:,1]-rho[:,0])/np.sqrt(nu)

    cnt = a.contourf(xv,tv,drho,50,cmap='RdBu')
    for c_ in cnt.collections:
        c_.set_edgecolor("face")

    a.set_xlabel(r'$x_\perp / L_{\rm Q}$')
    a.set_ylabel(r'$\omega_\varphi t$')
    a.set_xlim(-5.,5.)

    cb = f.colorbar(cnt,label=r'$\Pi_\varphi$',ticks = [-0.05,-0.025,0,0.025,0.05],pad=0.01)
    cb.solids.set_rasterized(True)
    
    return a.get_figure(), a

if __name__=="__main__":
    xv, tv, w, psi = load_data()
    tv_m = np.loadtxt('log.out',usecols=[4])
    nu = 0.01
    
    rho = np.sum(psi**2,axis=-2)
    drho = 0.5*(rho[:,1]-rho[:,0])/np.sqrt(nu)
    rhot = 0.5*(rho[:,0]+rho[:,1])/np.sqrt(nu)

    f,a = plt.subplots()
    cnt = a.contourf(xv,tv,rhot-rhot[0],50,cmap='RdBu', zorder=-10)
    a.set_rasterization_zorder(-1)

    #f,a = plot_rho_rel(xv,tv,psi,nu)
    #f,a = plot_rho_tot_deformation(xv,tv,psi)
    #a.text(0.95,0.95,r'$\bar{g} = 0.1$',va='top',ha='right',transform=a.transAxes)
    
    # This gets effective longitudinal field
    # Not debugged yet
    # rho = np.sum(psi**2,axis=-2)
    rho_l = np.sum(rho*w[np.newaxis,np.newaxis,:],axis=-1)
    drho_l = 0.5*(rho_l[:,1]-rho_l[:,0])/np.sqrt(nu) # Need to adjust w/ nu
    fk = np.fft.rfft(drho_l)

    #f,a = plt.subplots()
    #a.plot(tv,rho_l)
    
    #f,a = plt.subplots()
    #dom = 2.*np.pi/tv[-1]
    #om = dom*np.arange(fk.size)
    #a.plot(om,np.abs(fk)**2)
    #a.set_yscale('log')

    #f,a = plot_rho_tot_max(w,tv,psi)

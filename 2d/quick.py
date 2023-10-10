import numpy as np
import matplotlib.pyplot as plt

if __name__=="__main__":
    yv = np.loadtxt('ygrid.dat')
    xv = np.loadtxt('xgrid.dat')
    nu = 0.01
    
    ny = yv.size; nx = xv.size

    d = np.fromfile('fields.bin').reshape((-1,2,2,ny,nx))
    rho = np.sum(d**2,axis=2)
    drho = 0.5*(rho[:,1]-rho[:,0]) / np.sqrt(nu)
    rhot = np.mean(rho,axis=1) / np.sqrt(nu)

    dt = 0.625; dt_m = 2.*np.sqrt(nu)*dt
    xv_m = 2.*np.sqrt(nu)*xv
    yv_m = 2.*np.sqrt(nu)*yv

    drho0 = np.mean(drho,axis=(-2,-1))
    levs = np.linspace(-np.pi,np.pi,50)

    f,a = plt.subplots()
    for i_ in np.arange(drho.shape[0]):
        a.contourf(xv_m,yv_m,drho[i_]-drho0[i_],levs,extend='both',cmap='RdBu',zorder=-10)
        a.set_xlabel(r'$m_\varphi x$')
        a.set_ylabel(r'$m_\varphi y$')
        a.text(0.95,0.95,r'$m_\varphi t = %.2f$' % (dt_m*i_),ha='right',va='top',transform=a.transAxes)
        f.savefig('preheat-slice-2d-%05i.png' % i_)
        a.cla()
    

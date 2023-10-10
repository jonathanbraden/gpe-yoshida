#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

def multipanel_plot_2d(fName='fields_preheat_5term.bin', nLat=256, nu=0.01, slices = [0,1500,1550,1600,1650,1700,1800,1900,2000]):
    gs_dict = { 'wspace' : 0.02,
                'hspace' : 0.02
               }
    fig,a = plt.subplots(nrows=3, ncols=3, gridspec_kw=gs_dict)

    # Move all this into a separate data processing function
    f = np.fromfile(fName).reshape((-1,2,2,nLat,nLat))
    psi1 = f[:,0,0]+1j*f[:,0,1]
    psi2 = f[:,1,0]+1j*f[:,1,1]

    rho    = np.sum(f**2,axis=2)
    momRel = 0.5*(rho[:,1]-rho[:,0])/np.sqrt(nu)
    momRel_mean = np.mean(momRel,axis=(-2,-1))
    momRel = momRel - momRel_mean[:,np.newaxis,np.newaxis]
    # add dx and dt info

    exc = np.max(np.abs(momRel))
    levels = np.linspace(-0.9*exc,0.9*exc,51)
    for a_,s_ in zip(a.flatten(),slices):
        a_.contourf(momRel[s_],levels,cmap='RdBu',extend='both')
        a_.text(0.95,0.95,r'$m_\varphi t = $',va='top',ha='right',transform=a_.transAxes)
        a_.label_outer()

    fig.supxlabel(r'$m_\varphi x$')
    fig.supylabel(r'$m_\varphi y$')
    return fig,a

if __name__=="__main__":
    pass

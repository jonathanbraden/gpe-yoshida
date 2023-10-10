import numpy as np
import matplotlib.pyplot as plt

# Fix this to work on a periodic grid
def bright_soliton(x,t,eta,kappa,lSize):
    om = 0.5*(kappa**2-eta**2)
    return (eta/np.cosh(eta*(x-kappa*t)))*np.exp(1j*(kappa*x-om*t))

def make_snapshots(dat,xv,dt):
    for i,d_ in enumerate(dat):
        plt.fill_between(xv,-np.sqrt(d_[0,0,:]**2+d_[0,1,:]**2),np.sqrt(d_[0,0,:]**2+d_[0,1,:]**2),color='gray',alpha=0.2)
        plt.plot(xv,d_[0,0,:],'r-',alpha=0.8)
        plt.plot(xv,d_[0,1,:],'b-',alpha=0.9)
        plt.xlabel(r'$x$')
        plt.ylabel(r'$\psi$')
        plt.savefig('soliton-t%05i.png' %i)
        plt.clf()

def read_data(fName,n):
    return np.fromfile(fName).reshape((-1,2,2,n))
        
if __name__=="__main__":
    dat = read_data('fields.bin',1024)
    dx = 64./1024.
    dt = dx**2/4.*100
    make_snapshots(dat,np.linspace(0.,64.,1024,endpoint=False),dt)
                                                

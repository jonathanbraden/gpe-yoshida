import numpy as np
import matplotlib.pyplot as plt

# Can probably delete this
from scipy.linalg import circulant

# Move these derivative codes into a separate file
def laplacian_spectral_1d(f, dx):
    n = f.size
    norm = (2.*np.pi / dx)
    fk = -norm**2*np.fft.fft(f)

    # Check ordering
    df = fk * np.fft.fftfreq(n)**2

    return np.fft.ifft(df)

def grad_1d(f, dx):
    n = f.size
    norm = 2.*np.pi / dx
    fk = 1j*norm*np.fft.fft(f)

    df = fk*np.fft.fftfreq(n)

    return np.fft.ifft(df)

# Debug this
# Fix for non-square shapes
def laplacian_spectral(f, dx):
    nx, ny = f.shape[-2:]
    norm = (2.*np.pi / dx)
    fk = np.fft.fft2(f)

    # Check ordering
    df = -fk * (norm*np.fft.fftfreq(nx)[:,np.newaxis])**2
    df = df - fk * (norm*np.fft.fftfreq(ny)[np.newaxis,:])**2

    return np.fft.ifft2(df)

# Debug this
# Fix for non-square shapes
def grad_spectral(f, dx):
    nx, ny = f.shape[-2:]
    norm = 2.*np.pi / dx
    fk = 1j*norm*np.fft.fft2(f)

    # Check these orderings
    df_x = fk*np.fft.fftfreq(nx)[:,np.newaxis]
    df_y = fk*np.fft.fftfreq(ny)[np.newaxis,:]

    df_x = np.fft.ifft2(df_x)
    df_y = np.fft.ifft2(df_y)

    return df_x, df_y  # Merge into a single array


# Have to add correct normalization in here still
class Wavefunction():

    # Need to upgrade for non-square lattices
    def __init__(self, wf_file, pot_file, x_file, y_file):
        self.xv = np.loadtxt(x_file)
        self.yv = np.loadtxt(y_file)
        self.nx, self.ny = self.xv.size, self.yv.size

        self.nf = 1
        self.psi = np.fromfile(wf_file).reshape((-1,2*self.nf,self.ny,self.nx))
        self.psi = self.psi[:,0] + 1j*self.psi[:,1]

        self.pot = np.fromfile(pot_file).reshape((self.ny,self.nx))

        self.px = np.sum(np.abs(self.psi)**2,axis=-2)
        self.py = np.sum(np.abs(self.psi)**2,axis=-1)

        self.dx, self.dy = self.xv[1]-self.xv[0], self.yv[1]-self.yv[0]
        

    def get_prob(self):
        return np.abs(self.psi)**2*self.dx*self.dy
    
    def get_total_prob(self):
        return np.sum(np.abs(self.psi)**2*self.dx*self.dy,axis=(-2,-1))

    def marginal_densities(self, *, axis='x'):
        assert axis in ['x','y'], "Axis must be 'x' or 'y'"
        if (axis=='y'):
            return np.sum(np.abs(self.psi)**2, axis=-2)
        return np.sum(np.abs(self.psi)**2, axis=-1)

    def probability_current(self):
        return
    
    def reduced_density(self, tInd, axis):
        if (axis=='y'):
            pass
        elif (axis=='x'):
            pass
        else:
            print('Invalid axis')
            return None

    def von_neumann_entropy(self, axis='x'):
        return

    def shannon_entropy(self, axis):
        r"""
        Compute the Shannon entropy for the reduced probability distribution along the indicated axis.

        axis must be one of 'x' or 'y'
        """
        assert axis in ['x','y'], 'Error, axis must be x or y'

        if (axis=='x'):
            p = self.px / np.sum(self.px, axis=-1, keepdims=True)
        elif (axis=='y'):
            p = self.py / np.sum(self.py, axis=-1, keepdims=True)
        else:
            return None

        return -np.sum(p*np.log(p),axis=-1)

    # Write this, fix normalizations, etc.
    # Need to add some checks to avoid dividing by 0
    # Debug that I've added the mask correctly
    def kl_to_gaussian(self, axis='x',*, eps_g=1.e-7, eps_p=1.e-7, gauss_ref=False):
        """
        Compute KL Divergence from a Gaussian probability with same mean and variance
        """
        assert axis in ['x','y'], 'Error, axis must be x or y'
        
        if axis=='x':
            x,p = self.xv, self.px/np.sum(self.px,axis=-1,keepdims=True)
        elif axis=='y':
            x,p = self.yv, self.py/np.sum(self.py,axis=-1,keepdims=True)
        else:
            print(f'Error, axis must be x or y')
            return None

        mu = np.sum(x*p,axis=-1,keepdims=True)
        var = np.sum((x-mu)**2*p,axis=-1,keepdims=True)
        gauss = np.exp(-0.5*(x-mu)**2/var) / np.sqrt(2.*np.pi*var)
        gauss = gauss / np.sum(gauss,axis=-1,keepdims=True)
        
        mask = (np.abs(gauss) > eps_g) & (np.abs(p)>eps_p)

        if gauss_ref:
            kld = np.sum(gauss*np.log(gauss/p), axis=-1, where=mask)
        else:
            kld = np.sum(p*np.log(p/gauss), axis=-1, where=mask)
        
        return kld

    # Write this.  Makes minimal sense if the mean is oscillating
    # so I'll have to "shift the mean" somehow.
    def kl_from_initial(self):
        return
    
    def plot_potential(self, *, ax=None):
        r"""
        Make a contour plot of the 2-particle potential.
        """
        if ax is None:
            _, ax = plt.subplots()

        ax.contourf(self.xv, self.yv, self.pot)  # Check ordering
        return ax.get_figure(), ax

def compute_reduced_density_matrix(psi, method='matrix'):
    if method=='matrix':
        rho = _reduced_density_matrix(psi)
    elif method=='loop':
        rho = trace_density_matrix(psi)
    else:
        print(f'Error, invalid method')
        rho = None
    return rho
    
# This doesn't have dx norm part in it yet
# Improve speed
# Allow for a choice of which degree of freedom to trace out
def trace_density_matrix(psi):
    ny, nx = psi.shape  # Check if I want first or last index

    # Fix this
    # Oops, tracing out the zero mode
    rho = np.zeros((nx,nx), dtype=np.complex128)
    for i in range(nx):
        for j in range(nx):
            rho[i,j] = np.sum( np.conj(psi[i])*psi[j] ) # Fix this
    
    return rho

def reduced_density_matrix(psi, axis=0):
    if axis==0:
        rho_reduced = psi.T@np.conj(psi)
    else:
        rho_reduced = psi@np.conj(psi.T)

    return rho_reduced

from scipy.linalg import logm
def von_neumann_entropy(dens, *, method='eigen'):
    # Have to make sure density is normalized
    if (method=='eigen'):
        ev = np.sort(np.linalg.eigvalsh(dens))[::-1]
        entropy = -np.sum( ev[:20]*np.log(np.abs(ev[:20]) )) # check this
    if (method=='matrix'):
        entropy = -np.trace( dens*logm(dens) )
    return entropy

# Write code to compute Tr(rho^n) for replica trick
# Need matrix_pow command
def matrix_traces(dens, order):
    return
    
# Check the ordering of the directions in here
# Fix normalization
def probability_current(psi, dx):
    nt = psi.shape[0]
    nx, ny = psi.shape[-2:]
    current = np.zeros( psi.shape+(2,) )

    df_x, df_y = grad_spectral(psi, dx)
    current[...,0] = np.imag(np.conj(psi)*df_x)
    current[...,1] = np.imag(np.conj(psi)*df_y)
    
    return current

# Fix normalization here
def probability_divergence(psi, dx):
    lap = laplacian_spectral(psi, dx)
    
    return -np.imag(np.conj(psi)*lap)

# Write this to project transverse part of probability
def probability_transverse(psi, dx):
    return

# Write this to get scalar part of probability current
def probability_scalar(psi, dx):
    # Compute divergence (see above)
    # Then hit with inverse laplacian (but don't do the k=0 mode)

    # Then to get vector back, hit with derivative
    return

def movie_slides(prob,x,y,pot):
    fig, ax = plt.subplots()
    for i,p in enumerate(prob):
        ax.contour(x,y,pot,cmap='Reds')
        ax.contourf(x,y,p/np.max(p),levels=np.linspace(0.,1.,51),
                    extend='max',
                    cmap='binary',
                    alpha=0.75)
        fig.savefig('prob-%05i.png' % i)
        ax.cla()
    return
        
# Compute the effective potential of the mean field
def effective_potential_1d(psi, xv, yv, pot):
    for x_ in xv:
        pot_split = pot(x_+yv) + pot(x_-yv)

    # Now sum \int (V(x+r) + V(x-r)) \psi^2 dr
    # Figure out how to normalize out the psi(x) part
    return
        
if __name__=="__main__":
    wf = Wavefunction('fields.bin', 'trap.bin', 'xgrid.dat', 'ygrid.dat')
    
    prob = wf.get_prob()
    # This ordering is for the Fortran mapped to C ordering
    px = np.sum(prob,axis=-2)
    py = np.sum(prob,axis=-1)

    dprob = probability_divergence(wf.psi, wf.dx)
    current = probability_current(wf.psi, wf.dx)

    # Normalizing to density
    psi = wf.psi * np.sqrt(wf.dx*wf.dy)  # squared since two particles
    
    #entropy_mean = []
    #entropy_k = []
    #for p_ in psi:
    #    dens = p_.T@np.conj(p_)
    #    s = von_neumann_entropy(dens)
    #    entropy_mean.append(s)
    #    dens = p_@np.conj(p_).T
    #    s = von_neumann_entropy(dens)
    #    entropy_k.append(s)

    #eps = 1.e-5
    # I think the broadcasting on these is wrong
    #cond_prob_y = np.where(px>eps, prob/px, 0.) # Check broadcasting
    #cond_prob_x = np.where(py>eps, prob/py, 0.) # Check broadcasting
    #plt.plot(cond_prob_y[0])
    #plt.plot(py[0])

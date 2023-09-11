import numpy as np
import matplotlib.pyplot as plt

from scipy.linalg import circulant

# Have to add correct normalization in here still
class Wavefunction():

    # Need to upgrade for non-square lattices
    def __init__(self, wf_file, pot_file, x_file, y_file):
        self.xv = np.loadtxt(x_file)
        self.yv = np.loadtxt(y_file)
        nx, ny = self.xv.size, self.yv.size

        nf = 1
        self.psi = np.fromfile(wf_file).reshape((-1,2,ny,nx))
        self.psi = psi[:,0] + 1j*psi[:,1]

        self.pot = np.fromfield(pot_file).reshape((ny,nx))
        
        self.px = np.sum(np.abs(self.psi)**2,axis=-2)
        self.py = np.sum(np.abs(self.psi)**2,axis=-1)

    def density_matrix(self, tInd):
        return

    def marginal_densities(self, *, axis='x'):
        if (axis=='y'):
            return np.sum(np.abs(self.psi)**2, axis=-2)
        return np.sum(np.abs(self.psi)**2, axis=-1)

    def reduced_density(self, tInd, axis):
        if (axis=='y'):
            pass
        elif (axis=='x'):
            pass
        else:
            print('Invalid axis')
            return None
        
# This doesn't have dx norm part in it yet
# Improve speed
# Allow for a choice of which degree of freedom to trace out
def trace_density_matrix(psi):
    nx, ny = psi.shape  # Check if I want first or last index

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

def plot_re_im(psi,dt=0.25):
    nx,ny = psi.shape[-2:]
    tv = 0.25*np.arange(psi.shape[0])
    plt.plot(tv,psi[:,0,nx//2,ny//2])
    plt.plot(tv,psi[:,1,nx//2,ny//2])

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
    xv = np.loadtxt('xgrid.dat')
    yv = np.loadtxt('ygrid.dat')

    dx = xv[1]-xv[0]
    
    nx, ny = xv.size, yv.size
    n = nx  # Hack for now
    psi = np.fromfile('fields.bin').reshape((-1,2,n,n))
    psi = psi[:,0] + 1j*psi[:,1]

    pot = np.fromfile('trap.bin').reshape((n,n))
    prob = np.abs(psi)**2*dx**2
    # This ordering is for the Fortran mapped to C ordering
    px = np.sum(prob,axis=-2)
    py = np.sum(prob,axis=-1)

    #plt.contourf(xv,yv,pot)

    dprob = probability_divergence(psi, dx)
    current = probability_current(psi, dx)

    # Normalizing to density
    psi = psi * dx  # squared since two particles
    
    entropy_mean = []
    entropy_k = []
    for p_ in psi:
        dens = p_.T@np.conj(p_)
        s = von_neumann_entropy(dens)
        entropy_mean.append(s)
        dens = p_@np.conj(p_).T
        s = von_neumann_entropy(dens)
        entropy_k.append(s)

    eps = 1.e-5
    # I think the broadcasting on these is wrong
    #cond_prob_y = np.where(px>eps, prob/px, 0.) # Check broadcasting
    #cond_prob_x = np.where(py>eps, prob/py, 0.) # Check broadcasting
    #plt.plot(cond_prob_y[0])
    #plt.plot(py[0])

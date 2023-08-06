import numpy as np
import matplotlib.pyplot as plt

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
        ax.contourf(x,y,p,levels=np.linspace(0.,0.15,51),
                    extend='max',
                    cmap='binary',
                    alpha=0.75)
        fig.savefig('prob-%05i.png' % i)
        ax.cla()

if __name__=="__main__":
    xv = np.loadtxt('xgrid.dat')
    yv = np.loadtxt('ygrid.dat')

    dx = xv[1]-xv[0]
    
    n = xv.size
    psi = np.fromfile('fields.bin').reshape((-1,2,n,n))
    psi = psi[:,0] + 1j*psi[:,1]
    
    pot = np.fromfile('trap.bin').reshape((n,n))
    prob = np.abs(psi)**2

    #plt.contourf(xv,yv,prob[100]-prob[0])
    #plt.colorbar()

    #tv = 0.25*np.arange(101)
    #plt.plot(tv,psi[:,0,n//2,n//2])
    #plt.plot(tv,psi[:,1,n//2,n//2])

    #plt.contourf(xv,yv,pot)

    dprob = probability_divergence(psi, dx)
    current = probability_current(psi, dx)

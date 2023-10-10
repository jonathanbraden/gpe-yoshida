#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

def read_data(fName,nSims,nLat):
    return np.fromfile(fName).reshape((nSims,-1,nLat,nLat))

# Need to debug this (check odd lattices, etc)
# Then replace in various calls below
def wave_ind(n,full=True):
    """
    Convenience function to return wavenumber index along a single dimension.

    Input:
      n (integer) - Number of lattice sites
      full (Boolean) - If True then includes +ve and -ve wavenumbers
                       If False includes only +ve wavenumbers
    """
    nn = n//2
    if full:
        i_ind = np.concatenate( (np.arange(nn),np.arange(-nn,0)) )
    else:
        i_ind = np.arange(nn)
    return i_ind

# Assumes I've already ensemble averaged
# Bin this better (not floor binning)
def angle_average_spec_2d(fk,binFloor=True):
    nt = fk.shape[-3]    # Number of time-steps per sim
    nx = fk.shape[-2]    # Number of kx wavenumbers
    nnx = nx//2  # Check def of nnx for odd vs even
    nny = fk.shape[-1]   # Check def of nnx

    i_ind = wave_ind(nx,full=True)
    j_ind = np.arange(nny)
#    i_ind = np.concatenate( (np.arange(nnx),np.arange(-nnx,0)) )
# Check this is correct, then delete this commented line

    nk = np.floor(np.sqrt(nnx**2+nny**2)).astype(int) + 2
    power = np.zeros((nt,nk)); w = np.zeros(nk)

    if binFloor:
        for i_ in range(nx):
            ii_ = i_ind[i_]
            for j_ in range(nny):
                jj_ = j_ind[j_]
                k_ind = np.floor(np.sqrt(ii_**2+jj_**2)).astype(int)
                power[:,k_ind] += np.abs(fk[:,i_,j_])**2
                w[k_ind] += 1
    else:
        for i_ in range(nx):
            ii_ = i_ind[i_]
            for j_ in range(nny):
                jj_ = j_ind[j_]
                kv = np.sqrt(ii_**2+jj_**2)
                k_lt = np.floor(kv).astype(int); k_gt = k_lt+1
                # Vectorize the rest of this (test commented stuff
                #w_ = [w_lt, w_gt]  #insert expressions
                #power[:,k_lt:k_lt+2] += w_*np.abs(fk[:,i_,j_])**2
                #w[k_lt:k_lt+2] += w_
                w_lt = (1-(kv - k_lt)**2)**2
                w_gt = (1.-(k_gt-kv)**2)**2
                power[:,k_lt] += w_lt*np.abs(fk[:,i_,j_])**2
                power[:,k_gt] += w_gt*np.abs(fk[:,i_,j_])**2
                w[k_lt] += w_lt; w[k_gt] += w_gt
    #ii = np.where(count > 0)
    #power[:,ii] = power[:,ii]/count[ii]
    power = power / w
    return power

def angle_average_2d(fk,binFloor=False):
    """
    Given an input Fourier transform function, angle average it
    to produce the isotropized version.
    """
    nt = fk.shape[-3]    # Number of time-steps per sim
    nx = fk.shape[-2]    # Number of kx wavenumbers
    nnx = nx//2  # Check def of nnx for odd vs even
    nny = fk.shape[-1]   # Check def of nnx

    i_ind = wave_ind(nx,full=True)
    j_ind = np.arange(nny)
#    i_ind = np.concatenate( (np.arange(nnx),np.arange(-nnx,0)) )
# Check this is correct, then delete this commented line

    nk = np.floor(np.sqrt(nnx**2+nny**2)).astype(int) + 2
    fk_ave = np.zeros((nt,nk)); w = np.zeros(nk)

    if binFloor:
        for i_ in range(nx):
            ii_ = i_ind[i_]
            for j_ in range(nny):
                jj_ = j_ind[j_]
                k_ind = np.floor(np.sqrt(ii_**2+jj_**2)).astype(int)
                fk_ave[:,k_ind] += fk[:,i_,j_]
                w[k_ind] += 1
    else:
        for i_ in range(nx):
            ii_ = i_ind[i_]
            for j_ in range(nny):
                jj_ = j_ind[j_]
                kv = np.sqrt(ii_**2+jj_**2)
                k_lt = np.floor(kv).astype(int); k_gt = k_lt+1
                # Vectorize the rest of this (test commented stuff
                #w_ = [w_lt, w_gt]  #insert expressions
                #power[:,k_lt:k_lt+2] += w_*np.abs(fk[:,i_,j_])**2
                #w[k_lt:k_lt+2] += w_
                w_lt = (1-(kv - k_lt)**2)**2
                w_gt = (1.-(k_gt-kv)**2)**2
                fk_ave[:,k_lt] += w_lt*fk[:,i_,j_]
                fk_ave[:,k_gt] += w_gt*fk[:,i_,j_]
                w[k_lt] += w_lt; w[k_gt] += w_gt
    #ii = np.where(count > 0)
    #power[:,ii] = power[:,ii]/count[ii]
    fk_ave = fk_ave / w
    return fk_ave

# Input here is already a power spectrum
def angle_average_power_2d(pk):
    nt = pk.shape[-3]    # Number of time-steps per sim
    nx = pk.shape[-2]    # Number of kx wavenumbers
    nnx = nx//2  # Check def of nnx for odd vs even
    nny = pk.shape[-1]   # Check def of nnx
    i_ind = np.concatenate( (np.arange(nnx),np.arange(-nnx,0)) )
    j_ind = np.arange(nny)

    nk = np.floor(np.sqrt(nnx**2+nny**2)).astype(int) + 1
    power = np.zeros((nt,nk)); count = np.zeros(nk)
    for i_ in range(nx):
        ii_ = i_ind[i_]
        for j_ in range(nny):
            k_ind = np.floor(np.sqrt(ii_**2+j_**2)).astype(int)
            power[:,k_ind] += pk[:,i_,j_]
            count[k_ind] += 1
    #ii = np.where(count > 0)
    #power[:,ii] = power[:,ii]/count[ii]
    power = power / count
    return power

def angle_average_spec_1d(fk,kFull=True):
    """
    Angle average a 1D spectrum, assuming either the wavenumbers or the frequencies are purely positive.
    """
    return

def covariance_spectrum(phi,dphi):
    fk = np.fft.rfft2(phi)
    dfk = np.fft.rfft2(dphi)

    cov = [[np.mean(np.abs(fk)**2,axis=0), np.mean(fk*np.conj(dfk),axis=0)],
           [np.mean(fk*np.conj(dfk),axis=0), np.mean(np.abs(dfk)**2,axis=0)]]
    return np.array(cov)

def angle_average_cov_2d(cov):
    nt = cov.shape[-3]
    nx = cov.shape[-2]; nny = cov.shape[-1]
    nnx = nx//2+1
    i_ind = np.concatenate( (np.arange(nx//2),np.arange(-nx//2,0)) )
    j_ind = np.arange(nny) # Check this gives correct indices

    pow_size = np.floor(np.sqrt(nnx**2+nny**2)).astype(int) + 1
    power = np.zeros(nk); count = np.zeros(nk) # Fix this line
    for i_ in range(nx):
        ii_ = i_ind[i_]
        for j_ in range(nny):
            k_ind = np.floor(np.sqrt(ii_**2+j_**2)).astype(int) 
            power[:,k_ind] += cov[:,i_,j_]
            count[k_ind] += 1
    power = power / count  # do a where to avoid divide by zero
    return power

def get_k_vals(lSize,nk):
    dk = 2.*np.pi/lSize
    return dk*np.arange(nk)

def get_w_vals(tSize,nw):
    dw = 2.*np.pi/tSize
    return dw*np.arange(nw)

if __name__=="__main__":
    #inFile = 'data/quartic/len64_kc32_fixvar/l1-phi1/field-l1-phi1.bin'
    #inFile = 'data/double_well/field-dw-len64-kc32.bin'
    inFile = 'field.bin'; nS = 1; nL = 512
    f = np.fromfile(inFile).reshape((nS,-1,nL,nL))
    fk = np.fft.rfftn(f,axes=[-3,-2,-1])
    fk_t = np.fft.rfft2(f)
    # Fix this to not just use one realization
    pk = angle_average_spec_2d(fk[0],False)  # Welch filter
    pk_ = angle_average_2d(np.mean(np.abs(fk)**2,axis=0),False)
    
    pk_t = angle_average_spec_2d(fk_t[0],False)
    
    l = 64.; tf = 128.
    dk = 2.*np.pi/l; dw = 2.*np.pi/tf
    kv = dk*np.arange(pk.shape[-1])
    wv = dw*np.arange(pk.shape[0])

    # This is dependent on the particular model
    sig = 0.42
    dm2 = -0.5*sig*4.
    dm_g = (np.exp(-2.*sig)-1.)
    #lam = 1.
    #sig = 4.*0.0115*np.sqrt(lam)
    #dm2 = 3.*lam*sig
    #sig = 4.6e-2  #0.189*(1.5/1.)**2
    #dm2 = 0.5*sig*4.
    #dm4 = 3.*sig**2*16./24
    #dm6 = (3.*5.)*sig**3*2.**6/(24.*5.*6.)

    f,a = plt.subplots()

    import matplotlib as mpl
    import copy
    myCM = copy.copy(mpl.cm.get_cmap("OrRd"))
    myCM.set_under(alpha=0.)
    cnt = a.contourf(kv[:32],wv[:64],pk[:64,:32],21,vmin=0.5e9,extend='max',cmap=myCM)
    #myCM = copy.copy(mpl.cm.get_cmap("Blues"))
    #myCM.set_under(alpha=0.)
    #cnt = a.contourf(kv[:32],2.*wv[:64],pk[:64,:32],21,vmin=0.5e9,extend='max',cmap=myCM)
    a.set_ylim(0.,3.)
    a.plot(kv[:32],np.sqrt(kv[:32]**2+1.+dm2),'k--',alpha=0.5)
    a.plot(kv[:32],np.sqrt(kv[:32]**2+1.+dm_g),'r--',alpha=0.5)
    #plt.plot(kv[:32],np.sqrt(kv[:32]**2+1.-dm2+dm4),'r--')
    #plt.plot(kv[:32],np.sqrt(kv[:32]**2+1.-dm2+dm4-dm6),'g--')
    pass

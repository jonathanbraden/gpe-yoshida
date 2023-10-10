#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

class BEC_Data(object):
    def __init__(self,fName,nx,ny,nf):
        self.data = readData(fName,nx,ny,nf)
        self.nFld = nf

        self.psi = np.array([ self.data[:,i,0] + 1j*self.data[:,i,1] for i in range(nf) ])
        self.dphi = np.angle(np.conj(self.psi[0])*self.psi[1])
        self.theta = np.angle(self.psi[0]*self.psi[1])
        
def readData(fName,nx,ny,nf):
    return np.fromfile(fName).reshape((-1,nf,2,ny,nx))

def extract_relative_phase(psi2,psi1):
    x_sz = np.real(psi1*np.conj(psi2))
    y_sz = np.imag(np.conj(psi1)*psi2)
    return np.arctan2(y_sz,x_sz)

def gradient_squared(psi, dx, dy):
    nx = psi.shape[-1]; ny = psi.shape[-2]
    fk = np.fft.fft2(psi)
    
    kv = np.fft.fftfreq(nx); norm = 2.*np.pi/dx
    dpsi2 = np.abs( np.fft.ifft2(1j*kv[np.newaxis,:]*fk) )**2
    kv = np.fft.fftfreq(ny); norm = 2.*np.pi/dy
    dpsi2 = dpsi2 + np.abs( np.fft.ifft2(1j*kv[:,np.newaxis]*fk) )**2
    
    return dpsi2

if __name__=="__main__":
    pass

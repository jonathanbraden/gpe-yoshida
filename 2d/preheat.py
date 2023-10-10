import numpy as np
import matplotlib.pyplot as plt
from spectrum import angle_average_spec_2d

fName = 'fields_preheat_5term.bin'; nLat = 256
nu = 0.01
#fName = 'fields.bin'; nLat = 128

f = np.fromfile(fName).reshape((-1,2,2,nLat,nLat))
psi1 = f[:,0,0]+1j*f[:,0,1]
psi2 = f[:,1,0]+1j*f[:,1,1]

dphi = np.angle(np.conj(psi2)*psi1)
rho = np.sum(f**2,axis=2)
drho = 0.5*(rho[:,1]-rho[:,0])
rho_tot = 0.5*(rho[:,0]+rho[:,1])
theta = np.angle(psi1*psi2)

# Unwrap phi somewhere, or use cos?
dphi_mean = np.mean(dphi,axis=(-2,-1))

fk = np.fft.rfft2(drho)
pk = angle_average_spec_2d(fk)

#for i,p_ in enumerate(dphi):
#    plt.contourf(p_-dphi_mean[i],np.linspace(-np.pi,np.pi,50,endpoint=True),cmap='RdBu')
#    plt.savefig('preheat-no-mean-%05i.png' %i)
#    plt.clf()

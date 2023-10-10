import numpy as np
import matplotlib.pyplot as plt

if __name__=="__main__":
    d = np.fromfile('gpe_vacua.bin').reshape((5,2,-1))
    x = np.loadtxt('trap.dat',usecols=[0])

    rho = np.sum(d**2,axis=1)

    f,a = plt.subplots()
    for i_,r_ in enumerate(rho):
        a.plot(x,r_,label=r'$\bar{g} = 10^{%i}$' %  (i_-2))

    a.set_xlim(-6.5,6.5)
    a.set_ylabel(r'$\left|\psi_\perp\right|^2$')
    a.set_xlabel(r'$x_\perp / L_{\rm Q}$')
    a.legend()
    
    f.show()

#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

def load_data(fName,nlat,nf):
    return np.fromfile(fName).reshape((-1,nf,2,nlat))

def get_part_num(dat):
    return np.mean(dat[:,0,0]**2+dat[:,0,1]**2,axis=-1)

def part_num_convergence(files):
    f,a = plt.subplots()

    for f_ in files:
        dat = load_data(f_,16,2)
        rho = get_part_num(dat)
        a.plot(np.abs(rho-rho[0]))
    return f,a

def field_convergence(files):
    f,a = plt.subplots()
    nlat = 16; nfld = 2
    
    dat = [ load_data(f_,nlat,nfld) for f_ in files ]
    print(len(dat))
    for i_ in range(len(dat[:-1])):
        print(i_)
        print(dat[i_].shape)
        print(dat[i_+1].shape)
        a.plot(np.max(np.abs(dat[i_]-dat[i_+1]), axis=(1,2,3)))

    a.set_yscale('log')
    a.set_ylabel(r'$|\psi_{dt/2}-\psi_{dt}|_{\rm max}$')
    return f,a

def make_num_part_plots():
    f,a = part_num_convergence(files_o6)
    a.set_ylabel(r'$\rho-\rho_{\rm init}$')
    a.set_xlabel(r'time step')
    a.set_yscale('log')
    f.savefig('part-num-o6.pdf')

    f,a = part_num_convergence(files_o4)
    a.set_ylabel(r'$\rho-\rho_{\rm init}$')
    a.set_xlabel(r'time step')
    a.set_yscale('log')
    f.savefig('part-num-o4.pdf')

    f,a = part_num_convergence(files_o2)
    a.set_ylabel(r'$\rho-\rho_{\rm init}$')
    a.set_xlabel(r'time step')
    a.set_yscale('log')
    f.savefig('part-num-o2.pdf')
    return
        
if __name__=="__main__":
    files_o2 = ['fv-o2-dt0.01.bin','fv-o2-dt0.005.bin','fv-o2-dt0.0025.bin']
    files_o4 = ['fv-o4-dt0.01.bin','fv-o4-dt0.005.bin','fv-o4-dt0.0025.bin','fv-o4-dt0.00125.bin']
    files_o6 = ['fv-o6-dt0.01.bin','fv-o6-dt0.005.bin','fv-o6-dt0.0025.bin','fv-o6-dt0.00125.bin']

    f,a = field_convergence(files_o6)
    f.savefig('max-norm-o6.pdf')
    f,a = field_convergence(files_o4)
    f.savefig('max-norm-o4.pdf')
    f,a = field_convergence(files_o2)
    f.savefig('max-norm-o2.pdf')

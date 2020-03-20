import numpy as np
from matplotlib import pyplot as plt

def Sn(n,z):
    '''Calculates the following function:

    Sn(z,n) = z^n*sum_{Q\in Z}Q^{2n} exp(-z/2*Q^2) / sum_{Q\in Z}exp(-z/2*Q^2) 

    '''
    Qmax = 100
    num = 0.0
    denom = 1.0
    for Q in range(1,Qmax):
        num += 2*Q**(2*n)*np.exp(-0.5*z*Q**2)
        denom += 2*np.exp(-0.5*z*Q**2)
        
    return z**n*num/denom

if (__name__ == '__main__'):
    I = 0.25
    Tfinal = 4.0
    z = 4*np.pi**2*I/Tfinal
    eps=1.E-3
    zmax=12
    Z = np.arange(1.E-2,zmax,eps)
    S = np.vectorize(lambda x: Sn(1,x))
    R = np.vectorize(lambda x: 0.5*(Sn(2,x)-Sn(1,x)**2))
    do_plot = True
    if do_plot:
        plt.clf()
        ax = plt.gca()
        ax.set_xlabel('z')
        ax.set_xlim(0,zmax)
        ymin, ymax = 0, 1.75
        ax.set_ylim(ymin,ymax)
        plt.plot(Z,S(Z),
                 linewidth=2,
                 color='blue',
                 label='S(z)')
        plt.plot(Z,R(Z),
                 linewidth=2,
                 color='red',
                 label='R(z)')
        plt.plot([z,z],[ymin,ymax],
                 color='black',
                 linestyle='--',
                 linewidth=2)
        ax.legend(loc='upper right')
        plt.savefig('SandR.pdf',bbox_inches='tight')
    chi_t = 1/(4*np.pi**2*I)*S(z)
    Var_chi_t = 1/(8*np.pi**4*I**2)*R(z)
    print ('z          = '+('%10.6f' % z))
    print ('E[chi_t]   = '+('%10.6f' % chi_t))
    print ('Var[chi_t] = '+('%10.6f' % Var_chi_t))

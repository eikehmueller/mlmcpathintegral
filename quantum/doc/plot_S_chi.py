import numpy as np
from matplotlib import pyplot as plt

def S(z):
    nmax = 1000
    S1 = 0.
    S2 = 1.
    for Q in range(1,nmax):
        S1 += 2*Q**2*np.exp(-0.5*z*Q**2)
        S2 += 2*np.exp(-0.5*z*Q**2)
    return z*S1/S2

if (__name__ == '__main__'):
    T = 20.0
    I = 0.25
    z = 4*np.pi**2*I/T
    print (z, S(z), 1/(4.*np.pi**2*I)*S(z))
    eps = 1.E-2
    X = np.arange(eps,20,eps)
    Y = np.vectorize(S)(X)
    plt.clf()
    ax = plt.gca()
    ax.set_xlabel('z')
    ax.set_ylabel('S(z)')
    ax.set_ylim(0,1.1)
    plt.plot(X,Y,
             linewidth=2,
             color='blue')
    plt.savefig('S.pdf',bbox_inches='tight')

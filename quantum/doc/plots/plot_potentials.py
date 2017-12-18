import numpy as np
from matplotlib import pyplot as plt

def plot_potential(V,x_lim,y_lim,filename):
    plt.clf()
    ax = plt.gca()
    ax.set_xticks(())
    ax.set_yticks(())
    X = np.arange(x_lim[0],x_lim[1],0.01)
    Y = np.vectorize(V)(X)
    Y0 = 0.0*X-10.
    ax.set_xlim(x_lim[0],x_lim[1])
    ax.set_ylim(y_lim[0],y_lim[1])
    plt.axis('off')
    plt.plot(X,Y,linewidth=4,color='blue')
    plt.fill_between(X,Y,Y0,color='lightgray')
    plt.plot(x_lim,[0,0],linewidth=2,color='black',linestyle='--')
    plt.plot([0,0],y_lim,linewidth=2,color='black',linestyle='--')
    plt.savefig(filename,bbox_inches='tight')

def V_harmonic(x):
    return 0.5*x**2

def V_quartic(x):
    return -0.5*x**2+0.25*x**4

# Harmonic oscillator
x_lim = [-2,2]
y_lim = [-0.25,2]
filename = 'potential_harmonic.pdf'
plot_potential(V_harmonic,x_lim,y_lim,filename)

# Quartic oscillator
x_lim = [-2,2]
y_lim = [-0.5,2]
filename = 'potential_quartic.pdf'
plot_potential(V_quartic,x_lim,y_lim,filename)

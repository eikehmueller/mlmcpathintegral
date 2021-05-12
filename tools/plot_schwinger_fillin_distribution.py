import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

'''
Plot distribution of theta_j - theta_k for 4-dimensional
fillin distribution of quenched Schwinger model.

This is done both for the full distribution and for the Gaussian
approximation.
'''

def plot(ax,j,k,dataj,dataj_gaussian,datak,datak_gaussian):
    '''Plot histogram of theta_j - theta_k on given axis

    :arg ax: axis to plot on
    :arg j: index j
    :arg k: index k
    :dataj: dataframe with theta_j (full distribution)
    :dataj_gaussian: dataframe with theta_j (Gaussian approximation)
    :datak: dataframe with theta_k (full distribution)
    :datak_gaussian: dataframe with theta_k (Gaussian approximation)
    '''
    nbin = 64
    ax.hist(np.array(dataj)-np.array(datak),
            bins=nbin,
            density=True,
            alpha=0.5,
            label='full')
    ax.hist(np.array(dataj_gaussian)-np.array(datak_gaussian),
            bins=nbin,
            density=True,
            alpha=0.5,
            label='gaussian')
    ax.legend(loc='upper right')
    ax.set_title(r"$\theta_"+str(j)+r"-\theta_"+str(k)+"$",
                 y=0.8, pad=14)
    ax.set_xticks((-2*np.pi,-1.5*np.pi,-np.pi,-0.5*np.pi,0,+0.5*np.pi,+np.pi,+1.5*np.pi,+2*np.pi))
    ax.set_xticklabels(("$-2\pi$",r"$-\frac{3}{2}\pi$","$-\pi$",r"$-\frac{1}{2}\pi$","$0$",r"$+\frac{1}{2}\pi$","$+\pi$",r"$+\frac{3}{2}\pi$","$+2\pi$"))

#####################################################################
# M A I N
#####################################################################
if (__name__ == '__main__'):    
    filename = "../fillin_distribution.txt"
    filename_gaussian = "../fillin_distribution_gaussian.txt"
    data = pd.read_csv(filename,skiprows=5)
    data_gaussian = pd.read_csv(filename_gaussian,skiprows=5)

    plt.clf()
    fig, ax = plt.subplots(3,2,figsize=(14,10))
    ax[0,0].set_ylim(0,1.0)
    ax[0,1].set_ylim(0,1.0)

    plot(ax[0,0],1,2,
         data.theta1,data_gaussian.theta1,
         data.theta2,data_gaussian.theta2)
    plot(ax[1,0],1,3,
         data.theta1,data_gaussian.theta1,
         data.theta3,data_gaussian.theta3)
    plot(ax[2,0],1,4,
         data.theta1,data_gaussian.theta1,
         data.theta4,data_gaussian.theta4)
    plot(ax[0,1],2,3,
         data.theta2,data_gaussian.theta2,
         data.theta3,data_gaussian.theta3)
    plot(ax[1,1],2,4,
         data.theta2,data_gaussian.theta2,
         data.theta4,data_gaussian.theta4)
    plot(ax[2,1],3,4,
         data.theta3,data_gaussian.theta3,
         data.theta4,data_gaussian.theta4)
    plt.savefig("schwinger_fillin_distribution.pdf",bbox_inches="tight")

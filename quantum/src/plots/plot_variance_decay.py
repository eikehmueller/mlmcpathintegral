'''Theoretical performance Analysis of MonAte Carlo methods'''
import sys
import subprocess
import math
import numpy as np
from multilevel import *

def plot_variance_decay(filenames,labels,output_filename,power):
    '''Plot variance of differences as a function of the lattice size
        '''
    plt.clf()
    p_var = []
    l_var = []
    markers=['o','s']
    marker_facecolors=['white','blue']
    i = 0
    for filename, label in zip(filenames,labels):
        M_lat = []
        variance = []
        diff_variance = []
        with open(filename) as f:
            for line in f.readlines():
                if (not (re.match(' *#.*',line))):
                    M_lat_tmp, variance_tmp, diff_variance_tmp = line.split()
                    M_lat.append(int(float(M_lat_tmp)))
                    variance.append(float(variance_tmp))
                    diff_variance.append(float(diff_variance_tmp))
        M_lat = np.array(M_lat)
        variance = np.array(variance)
        diff_variance = np.array(diff_variance)
        ax = plt.gca()
        ax.set_ylim((1.E-5,1000.))
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'number of time slices $M_\ell$')
        ax.set_ylabel(r'Variance')
        plt.plot(M_lat,diff_variance,
                 linewidth=0,
                 color='blue',
                 marker=markers[i],
                 markeredgewidth=2,
                 markeredgecolor='blue',
                 markerfacecolor=marker_facecolors[i])[0]
        p_var.append(plt.plot([1E9,1E9],[1E9,1E9],
                              linewidth=2,
                              color='blue',
                              marker=markers[i],
                              markeredgewidth=2,
                              markeredgecolor='blue',
                              markerfacecolor=marker_facecolors[i])[0])
        l_var.append(r'Var[$Y_\ell$] ('+label+')')
        p_var0 = plt.plot(M_lat,variance,
                          linewidth=2,
                          color='black',
                          linestyle='--',
                          marker='s',
                          markeredgewidth=2,
                          markeredgecolor='black',
                          markerfacecolor='white')[0]
        ax.set_xlim(0.5*M_lat[0],2.*M_lat[-1])
        ax.set_xticks(M_lat)
        ax.set_xticklabels([str(x) for x in M_lat])
        M_ref = M_lat[0]
        var_ref = diff_variance[0]
        rho = 8.
        # fit to data
        log_variance = np.log(diff_variance)
        log_M_lat = np.log(M_lat)
        a_1,a_0 = np.polyfit(log_M_lat,log_variance,1)
        M_lat = np.exp(np.arange(math.log(M_lat[0])-0.25,math.log(M_lat[-1])+0.25,0.01))
        C_deltaV = math.exp(a_0)
        p_fit = plt.plot(M_lat,C_deltaV*M_lat**a_1,
                         linewidth=2,
                         color='blue')[0]
        i += 1

        
    p_guide = plt.plot([M_ref,rho*M_ref],[0.5*var_ref,0.5/rho**power*var_ref],
                       linewidth=2,
                       color='green')[0]
    plt.legend(p_var+[p_var0,p_guide],
               l_var+[r'Var[$X^2$]',
                r'$\propto M_\ell^{-'+str(power)+'}$'],'upper left',ncol=2,fontsize=18)

    plt.savefig(output_filename,bbox_inches='tight')

    
if (__name__ == '__main__'):
    # *** harmonic oscillator ***

    T_final = 1.0
    m0 = 1.0
    mu2 = 1.0
    lmbda = 1.0

    action = 'harmonic oscillator'

    print " action          = ", action
    print " m0              = ", m0
    print " mu^2            = ", mu2
    print " lambda          = ", lmbda
    print " T_final         = ", T_final
    print

    do_run = False

    executable = 'driver_ho.x'
    filenames = []
    labels = []
    multilevel = MultiLevel(action,T_final,m0,mu2,lmbda,0,'driver_ho.x')
    filename='variance_ho_not_renormalised.dat'
    filenames.append(filename)
    labels.append('naive')
    if (do_run):
        multilevel._run_all(filename=filename,force=True)
    multilevel = MultiLevel(action,T_final,m0,mu2,lmbda,2,'driver_ho.x')
    filename = 'variance_ho_renormalised.dat'
    if (do_run):
        multilevel._run_all(filename=filename,force=True)
    filenames.append(filename)
    labels.append('renormalised')
    plot_variance_decay(filenames,labels,'variance_decay_ho.pdf',2)

    # *** quartic oscillator ***
    do_run = True

    T_final = 1.0
    m0 = 1.0
    mu2 = -1.0
    lmbda = 1.0

    action = 'quartic oscillator'

    print " action          = ", action
    print " m0              = ", m0
    print " mu^2            = ", mu2
    print " lambda          = ", lmbda
    print " T_final         = ", T_final
    print

    executable = 'driver_qo.x'
    filenames = []
    labels = []
    multilevel = MultiLevel(action,T_final,m0,mu2,lmbda,0,'driver_qo.x')
    filename='variance_qo_not_renormalised.dat'
    filenames.append(filename)
    labels.append('naive')
    if (do_run):
        multilevel._run_all(filename=filename,force=True)
    plot_variance_decay(filenames,labels,'variance_decay_qo.pdf',1)

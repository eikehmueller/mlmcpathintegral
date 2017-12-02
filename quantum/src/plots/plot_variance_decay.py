import sys
import re
import subprocess
import math
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({'font.size': 22})

executable = "../c/driver.x"

m0 = 1.0
mu2 = 1.0
T_final = 1.0
n_burnin = 10000
n_samples = 1000000

def run_variance(M_lat):
    '''Run code and return variance
    
    :arg M_lat: Number of lattice points
    '''
    parameters = '''
    M_lat           %(M_LAT)d
    T_final         %(T_FINAL)f
    m0              %(M0)f
    mu2             %(MU2)f
    n_burnin        %(N_BURNIN)d
    n_samples       %(N_SAMPLES)d
    ''' % {'M_LAT': M_lat,
           'T_FINAL':T_final,
           'M0':m0,
           'MU2':mu2,
           'N_BURNIN':n_burnin,
           'N_SAMPLES':n_samples}
    with open('parameters.in','w') as f:
        print >> f, parameters 
    output = subprocess.check_output([executable])
    for line in output.split('\n'):
        m = re.match('^ *mean *= *(.*) *variance *(.*)$',line)
        if m:
            mean = float(m.group(1))
            variance = float(m.group(2))
    return variance

def run_all(M_lat_list):
    with open('variance.dat','w') as f:
        for M_lat in M_lat_list:
            variance = run_variance(M_lat)
            print >> f, ('%12.8f' % M_lat)+'  '+('%12.8f' % variance)

def plot_all():
    M_lat = []
    variance = []
    with open('variance.dat') as f:
        for line in f.readlines():
            M_lat_tmp, variance_tmp = line.split()
            M_lat.append(int(float(M_lat_tmp)))
            variance.append(float(variance_tmp))
    M_lat = np.array(M_lat)
    variance = np.array(variance)
    plt.clf()
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'number of time slices $M_\ell$')
    ax.set_ylabel(r'Var($Y_\ell$)')
    p_var = plt.plot(M_lat,variance,
                     linewidth=2,
                     color='blue',
                     marker='o',
                     markeredgewidth=2,
                     markeredgecolor='blue',
                     markerfacecolor='white')[0]
    ax.set_xlim(0.5*M_lat[0],2.*M_lat[-1])
    ax.set_xticks(M_lat)
    ax.set_xticklabels([str(x) for x in M_lat])
    M_ref = M_lat[0]
    var_ref = variance[0]
    rho = 8.
    p_lin = plt.plot([M_ref,rho*M_ref],[0.5*var_ref,0.5/rho*var_ref],
                     linewidth=2,
                     color='black')[0]
    # Linear fit to data
    log_variance = np.log(variance)
    log_M_lat = np.log(M_lat)
    a_1,a_0 = np.polyfit(log_M_lat,log_variance,1)
    M_lat = np.exp(np.arange(math.log(M_lat[0])-0.25,math.log(M_lat[-1])+0.25,0.01))
    C_deltaV = math.exp(a_0)
    p_fit = plt.plot(M_lat,C_deltaV*M_lat**a_1,
                     linewidth=2,
                     color='red')[0]
    print 'C_deltaV = ',C_deltaV

    plt.legend((p_var,p_lin,p_fit),('variance','linear decay','fit'),'upper right')
    plt.savefig('variance_decay.pdf',bbox_inches='tight')

if (__name__ == '__main__'):
    do_run = False
    do_plot = True
    M_lat_list = (8,16,32,64,128)
    if (do_run):
        run_all(M_lat_list)
    if (do_plot):
        plot_all()

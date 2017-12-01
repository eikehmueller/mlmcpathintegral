import sys
import re
import subprocess
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
            a_lat = T_final/M_lat
            variance = run_variance(M_lat)
            print >> f, ('%12.8f' % a_lat)+'  '+('%12.8f' % variance)

def plot_all():
    a_lat = []
    variance = []
    with open('variance.dat') as f:
        for line in f.readlines():
            a_lat_tmp, variance_tmp = line.split()
            a_lat.append(float(a_lat_tmp))
            variance.append(float(variance_tmp))
    a_lat = np.array(a_lat)
    variance = np.array(variance)
    plt.clf()
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'lattice spacing $a_\ell$')
    ax.set_ylabel(r'Var($Y_\ell$)')
    p_var = plt.plot(a_lat,variance,
                     linewidth=2,
                     color='blue',
                     marker='o',
                     markeredgewidth=2,
                     markeredgecolor='blue',
                     markerfacecolor='white')[0]
    a_ref = a_lat[0]
    var_ref = variance[0]
    p_lin = plt.plot([0.2*a_ref,a_ref],[0.5*0.2*var_ref,0.5*var_ref],
                     linewidth=2,
                     color='black')[0]
    plt.legend((p_var,p_lin),('variance','linear decay'),'upper left')
    plt.savefig('variance_decay.pdf',bbox_inches='tight')

if (__name__ == '__main__'):
    do_run = False
    do_plot = True
    M_lat_list = (8,16,32,64,128)
    if (do_run):
        run_all(M_lat_list)
    if (do_plot):
        plot_all()

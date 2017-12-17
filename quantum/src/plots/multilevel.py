'''Cost parameters of multilevel method'''
import sys
import subprocess
import math
import numpy as np
import re
import os
from matplotlib import pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({'font.size': 18})

class MultiLevel(object):
    '''Class for analysing and plotting multilevel data

    :arg action: Name of action
    :arg T_final: Final time
    :arg m0: Particle mass m_0
    :arg mu2: Quadratic potential parameter \mu^2
    :arg lambda: Quartic potential parameter \lambda
    '''
    def __init__(self,action,T_final,m0,mu2,lmbda):
        self.action = action
        self.T_final = T_final
        self.m0 = m0
        self.mu2 = mu2
        self.lmbda = lmbda
        self.perturbative = 0
        self.n_burnin = 10000
        self.n_samples = 100000
        self.hmc_sampling = 0
        if (not (self.action=='harmonic oscillator')):
            self.hmc_sampling = 1
        self.T_hmc = 40.0
        self.dt_hmc = 0.01
        self.n_burnin_hmc  = 10000
        
        self.M_lat_list = (8,16,32,64,128)
        self._C_deltaV = None
        self.executable = "../c/driver.x"

    def _run_variance(self,M_lat):
        '''Run code and return variance
        
        :arg M_lat: Number of lattice points
        '''
        parameters = '''
        M_lat           %(M_LAT)d
        T_final         %(T_FINAL)f
        m0              %(M0)f
        mu2             %(MU2)f
        lambda          %(LAMBDA)f
        perturbative    %(PERTURBATIVE)d
        n_burnin        %(N_BURNIN)d
        n_samples       %(N_SAMPLES)d
        hmc_sampling    %(HMC_SAMPLING)d
        T_hmc           %(T_HMC)f
        dt_hmc          %(DT_HMC)f
        n_burnin_hmc    %(N_BURNIN_HMC)d
        ''' % {'M_LAT': M_lat,
               'T_FINAL':self.T_final,
               'M0':self.m0,
               'MU2':self.mu2,
               'LAMBDA':self.lmbda,
               'PERTURBATIVE':self.perturbative,
               'N_BURNIN':self.n_burnin,
               'N_SAMPLES':self.n_samples,
               'HMC_SAMPLING':self.hmc_sampling,
               'T_HMC':self.T_hmc,
               'DT_HMC':self.dt_hmc,
               'N_BURNIN_HMC':self.n_burnin_hmc}
        with open('parameters.in','w') as f:
            print >> f, parameters 
        output = subprocess.check_output([self.executable])
        for line in output.split('\n'):
            m = re.match('^ *Action *= *(.*) *$',line)
            if m:
                action=m.group(1).strip()
            m = re.match('^ *mean *= *(.*) *variance *(.*)$',line)
            if m:
                mean = float(m.group(1))
                diff_variance = float(m.group(2))
            m = re.match('^ *Var\[x\^2\] *= *(.*) *$',line)
            if m:
                variance = float(m.group(1))
        if (not (action==self.action)):
            print " ERROR: Action mismatch. Found \'",action,"\' expected \'",self.action+"\'"
        return variance, diff_variance

    def _run_all(self):
        '''
        Run code for a range of different lattice size and save both
        the lattice size and the variance to a file
        '''
        # Check if file already exists and was generated for the same
        # parameters
        tolerance = 1.E-12
        recreate_data = False
        match_T_final = False
        match_m0 = False
        match_mu2 = False
        recreate_data = not (os.path.isfile('variance.dat'))
        if (recreate_data):
            with open('variance.dat','w') as f:
                print >> f, "# T_final      = ",self.T_final
                print >> f, "# m0           = ",self.m0
                print >> f, "# mu2          = ",self.mu2
                print >> f, "# lambda       = ",self.lmbda
                print >> f, "# pertubative  = ",self.perturbative
                print >> f, "# n_burnin     = ",self.n_burnin
                print >> f, "# n_samples    = ",self.n_samples
                print >> f, "# hmc_sampling = ",self.hmc_sampling
                print >> f, "# T_hmc        = ",self.T_hmc
                print >> f, "# dt_hmc       = ",self.dt_hmc
                print >> f, "# n_burnin_hmc = ",self.n_burnin_hmc

                for M_lat in self.M_lat_list:
                    variance, diff_variance = self._run_variance(M_lat)
                    print >> f, ('%12.8d' % M_lat)+'  '+('%12.8f' % variance)+'  '+('%12.8f' % diff_variance)

    def plot_variance_decay(self):
        '''Plot variance of differences as a function of the lattice size
        '''
        self._run_all()
        M_lat = []
        variance = []
        diff_variance = []
        with open('variance.dat') as f:
            for line in f.readlines():
                if (not (re.match(' *#.*',line))):
                    M_lat_tmp, variance_tmp, diff_variance_tmp = line.split()
                    M_lat.append(int(float(M_lat_tmp)))
                    variance.append(float(variance_tmp))
                    diff_variance.append(float(diff_variance_tmp))
        M_lat = np.array(M_lat)
        variance = np.array(variance)
        diff_variance = np.array(diff_variance)
        plt.clf()
        ax = plt.gca()
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'number of time slices $M_\ell$')
        ax.set_ylabel(r'Var($Y_\ell$)')
        p_var = plt.plot(M_lat,diff_variance,
                         linewidth=2,
                         color='blue',
                         marker='o',
                         markeredgewidth=2,
                         markeredgecolor='blue',
                         markerfacecolor='white')[0]
        p_var0 = plt.plot(M_lat,variance,
                          linewidth=2,
                          color='black',
                          linestyle='--',
                          marker='o',
                          markeredgewidth=2,
                          markeredgecolor='black',
                          markerfacecolor='white')[0]
        ax.set_xlim(0.5*M_lat[0],2.*M_lat[-1])
        ax.set_xticks(M_lat)
        ax.set_xticklabels([str(x) for x in M_lat])
        M_ref = M_lat[0]
        var_ref = variance[0]
        rho = 8.
        p_lin = plt.plot([M_ref,rho*M_ref],[0.5*var_ref,0.5/rho*var_ref],
                         linewidth=2,
                         color='green')[0]
        # Linear fit to data
        log_variance = np.log(diff_variance)
        log_M_lat = np.log(M_lat)
        a_1,a_0 = np.polyfit(log_M_lat,log_variance,1)
        M_lat = np.exp(np.arange(math.log(M_lat[0])-0.25,math.log(M_lat[-1])+0.25,0.01))
        self._C_deltaV = math.exp(a_0)
        p_fit = plt.plot(M_lat,self._C_deltaV*M_lat**a_1,
                         linewidth=2,
                         color='red')[0]

        plt.legend((p_var,p_var0,p_fit,p_lin),
                   (r'Var[$Y_\ell$]',
                    r'Var[$\langle X^2\rangle$]',
                    'linear fit',
                    r'$\propto M_\ell^{-1}$'),'upper left',ncol=2,fontsize=18)
        plt.savefig('variance_decay.pdf',bbox_inches='tight')

    @property
    def C_deltaV(self):
        if (self._C_deltaV == None):
            self.plot_variance_decay()
        return self._C_deltaV

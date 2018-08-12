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
    :arg renormalisation: Renormalisation to use
    :arg executable: Name of executable
    '''
    def __init__(self,action,T_final,m0,mu2,lmbda,renormalisation,executable):
        self.action = action
        self.T_final = T_final
        self.m0 = m0
        self.mu2 = mu2
        self.lmbda = lmbda
        self.sigma = 1.0
        self.renormalisation = renormalisation
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
        self.executable = "../c/"+executable

    def _run_variance(self,M_lat,filename,force):
        '''Run code and return variance
        
        :arg M_lat: Number of lattice points
        :arg filename: Save to specific filename
        :arg force: Run in any case?
        '''
        parameters = '''
        M_lat           %(M_LAT)d
        T_final         %(T_FINAL)f
        m0              %(M0)f
        mu2             %(MU2)f
        lambda          %(LAMBDA)f
        sigma           %(SIGMA)f
        renormalisation %(RENORMALISATION)d
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
               'SIGMA':self.sigma,
               'RENORMALISATION':self.renormalisation,
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
            m = re.match('^ *Var\[QoI\] *= *(.*) *$',line)
            if m:
                variance = float(m.group(1))
        if (not (action==self.action)):
            print " ERROR: Action mismatch. Found \'",action,"\' expected \'",self.action+"\'"
        return variance, diff_variance

    def _run_all(self,filename='variance.dat',force=False):
        '''
        Run code for a range of different lattice size and save both
        the lattice size and the variance to a file

        :arg filename: Save to specific filename
        :arg force: Run in any case?
        '''
        # Check if file already exists and was generated for the same
        # parameters
        tolerance = 1.E-12
        recreate_data = False
        match_T_final = False
        match_m0 = False
        match_mu2 = False
        recreate_data = (not (os.path.isfile(filename))) or force
        if (recreate_data):
            with open(filename,'w') as f:
                print >> f, "# T_final         = ",self.T_final
                print >> f, "# m0              = ",self.m0
                print >> f, "# mu2             = ",self.mu2
                print >> f, "# lambda          = ",self.lmbda
                print >> f, "# sigma           = ",self.sigma
                print >> f, "# renormalisation = ",self.renormalisation
                print >> f, "# n_burnin        = ",self.n_burnin
                print >> f, "# n_samples       = ",self.n_samples
                print >> f, "# hmc_sampling    = ",self.hmc_sampling
                print >> f, "# T_hmc           = ",self.T_hmc
                print >> f, "# dt_hmc          = ",self.dt_hmc
                print >> f, "# n_burnin_hmc    = ",self.n_burnin_hmc

                for M_lat in self.M_lat_list:
                    variance, diff_variance = self._run_variance(M_lat,filename,force)
                    print >> f, ('%12.8d' % M_lat)+'  '+('%12.8f' % variance)+'  '+('%12.8f' % diff_variance)

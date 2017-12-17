'''Theoretical performance Analysis of Monte Carlo methods'''
import sys
import subprocess
import math
import numpy as np
import re
import os
from matplotlib import pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({'font.size': 18})
from singlelevel import *
from multilevel import *

class CostAnalysis(object):
    '''
    Class for comparative cost analysis

    :arg action: Name of action
    :arg T_final: Final time
    :arg m0: Particle mass m_0
    :arg mu2: Quadratic potential parameter \mu^2
    '''
    def __init__(self,action,T_final,m0,mu2):
        self.action = action
        self.T_final = T_final
        self.m0 = m0
        self.mu2 = mu2
        self.M0 = 8
        self.singlelevel = SingleLevel(self.T_final,self.m0,self.mu2)
        self.multilevel = MultiLevel(self.action,self.T_final,self.m0,self.mu2,0.0)

    def cost_singlelevel(self,epsilon):
        '''Theoretical cost estimate for single level method

        :arg epsilon: Tolerance on root mean square error
        '''
        C_bias = abs(self.singlelevel.Xsquared_bias_constant())
        variance = 2.*self.singlelevel.Xsquared_continuum()**2
        return 2.**(5./4.)*variance*math.sqrt(C_bias)*epsilon**(-2.5)

    def cost_multilevel(self,epsilon):
        '''Theoretical cost estimate for multi level method

        :arg epsilon: Tolerance on root mean square error
        '''
        C_bias = self.singlelevel.Xsquared_bias_constant()
        variance = 2.*self.singlelevel.Xsquared_continuum()**2
        C_deltaV = self.multilevel.C_deltaV
        L_level = 0.25*(1.+math.log(C_bias**2*self.M0**(-4)))/math.log(2.) - 0.5*math.log(epsilon)/math.log(2.)
        return 4.*C_deltaV*epsilon**(-2)*(L_level + math.sqrt(self.M0*variance/(2.*C_deltaV)))**2
        
    def plot_cost(self):
        '''Plot the cost as a function of the lattice spacing.

        If M is the number of time slices and N the number of MC samples, 
        then balancing the bias- and sampling-error gives a cost of
        '''
        M_lat = np.array([8,16,32,64,128,256,512])

        # Single level cost
        C_bias = abs(self.singlelevel.Xsquared_bias_constant())
        epsilon = math.sqrt(2.)*C_bias*(1.*M_lat)**(-2)
        cost_singlelevel = np.vectorize(self.cost_singlelevel)(epsilon)
        cost_multilevel = np.vectorize(self.cost_multilevel)(epsilon)
        
        plt.clf()
        ax = plt.gca()
        ax.set_xscale('log')
        ax.set_yscale('log')
        p_single = plt.plot(epsilon,cost_singlelevel,
                            linewidth=2,
                            color='blue',
                            marker='o',
                            markersize=8)[0]
        p_multi = plt.plot(epsilon,cost_multilevel,
                           linewidth=2,
                           color='red',
                           marker='o',
                           markersize=8)[0]
        eps_ref = epsilon[-1]
        cost_ref = cost_singlelevel[-1]
        rho = 10.0
        p_theory = plt.plot([eps_ref,rho*eps_ref],
                            [0.25*cost_ref,0.25*rho**(-2.5)*cost_ref],
                            linewidth=2,
                            color='black')[0]
        ax.set_xlabel(r'tolerance $\epsilon$')
        ax.set_ylabel(r'Cost')
        plt.legend((p_single,p_multi,p_theory),
                   (r'$\operatorname{Cost}^{\operatorname{StMC}}$',
                    r'$\operatorname{Cost}^{\operatorname{MLMC}}$',
                    r'$\propto \epsilon^{-2.5}$'),
                   'upper right')
        plt.savefig('cost.pdf',bbox_inches='tight')

    
if (__name__ == '__main__'):
    # Numerical values for parameters
    T_final = 1.0
    m0 = 1.0
    mu2 = 1.0
    lmbda = 1.0

    action = 'harmonic_oscillator'

    print " action  = ", action
    print " m0      = ", m0
    print " mu^2    = ", mu2
    print " lambda  = ", lmbda
    print " T_final = ", T_final
    print

    
    singlelevel = SingleLevel(T_final,m0,mu2)
    print " C_{bias} = ",singlelevel.Xsquared_bias_constant()
    print " V        = ",2.*singlelevel.Xsquared_continuum()**2
    singlelevel.plot_bias()
    
    multilevel = MultiLevel(action,T_final,m0,mu2,lmbda)
    multilevel.plot_variance_decay()
    print "delta C_V = ",multilevel.C_deltaV
    
    cost_analysis = CostAnalysis(action,T_final,m0,mu2)
    cost_analysis.plot_cost()


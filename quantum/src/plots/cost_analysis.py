'''Theoretical performance Analysis of Monte Carlo methods'''
import sys
import subprocess
import math
import numpy as np
import re
from matplotlib import pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({'font.size': 18})

class SingleLevel(object):
    def __init__(self,T_final,m0,mu2):
        self.T_final = T_final
        self.m0 = m0
        self.mu2 = mu2
        
    def Xsquared(self,M_lat):
        '''Calculate sigma^2 = <X^2> for a given number of
        timeslices.
        
        :arg M_lat: Number of timeslices
        '''
        a_lat = self.T_final/M_lat
        R = 1.+0.5*self.mu2*a_lat**2 - a_lat*math.sqrt(self.mu2)*math.sqrt(1.+0.25*a_lat**2*self.mu2)
        return 1./(2.*self.m0*math.sqrt(self.mu2)*math.sqrt(1+0.25*a_lat**2*self.mu2))*(1.+R**M_lat)/(1.-R**M_lat)

    def Xsquared_continuum(self):
        '''Calculate continuum limit of sigma^2 = <X^2>'''
        R = math.exp(-math.sqrt(self.mu2)*self.T_final)
        return 1./(2.*self.m0*math.sqrt(self.mu2))*(1.+R)/(1.-R)


    def Xsquared_bias_constant(self):
        ''' The bias of <X^2> is C_{bias}*M^{-2} 
        This method calculates C_{bias} by numerical differentiation
        '''
        M_ref = 1024
        numerical = (self.Xsquared(M_ref)-self.Xsquared_continuum())*M_ref**2
        tmp = 0.5*math.sqrt(self.mu2)*self.T_final
        analytical = math.sqrt(self.mu2)*self.T_final**2/(16.*self.m0)*(-1./3.*tmp/math.sinh(tmp)**2 + 1./math.tanh(tmp))
        return abs(numerical)
    
    def plot_bias(self):
        '''Produce plot of bias for a range of different grid sizes'''
        M_lat = np.array([8,16,32,64,128,256])
        bias = np.vectorize(lambda M: abs(self.Xsquared(M)-self.Xsquared_continuum()))(M_lat)
        plt.clf()
        ax = plt.gca()
        ax.set_xscale('log')
        ax.set_yscale('log')
        p_bias = plt.plot(M_lat,bias,
                          linewidth=2,
                          color='blue',
                          marker='o',
                          markersize=8)[0]
        M_ref = M_lat[0]
        bias_ref = bias[0]
        rho = 4.0
        p_quadratic = plt.plot([M_ref,rho*M_ref],
                               [0.5*bias_ref,0.5*rho**(-2)*bias_ref],
                               linewidth=2,
                               color='black')[0]
        ax.set_xlabel(r'number of time slices $M$')
        ax.set_ylabel(r'bias of $\langle X^2\rangle$')
        ax.set_xticks(M_lat)
        ax.set_xticklabels([str(x) for x in M_lat])
        plt.legend((p_bias,p_quadratic),
                   ('bias','quadratic decay'),
                   'upper right')
        plt.savefig('bias.pdf',bbox_inches='tight')

class MultiLevel(object):
    '''Class for analysing and plotting multilevel data

    :arg T_final: Final time
    :arg m0: Particle mass m_0
    :arg mu2: Quadratic potential parameter \mu^2
    '''
    def __init__(self,T_final,m0,mu2):
        self.T_final = T_final
        self.m0 = m0
        self.mu2 = mu2
        self.perturbative = 0
        self.n_burnin = 10000
        self.n_samples = 1000000
        self.M_lat_list = (8,16,32,64,128)
        self._C_deltaV = None
        self.executable = "../c/driver.x"
        self._singlelevel = SingleLevel(self.T_final,self.m0,self.mu2)

    def _run_variance(self,M_lat):
        '''Run code and return variance
        
        :arg M_lat: Number of lattice points
        '''
        parameters = '''
        M_lat           %(M_LAT)d
        T_final         %(T_FINAL)f
        m0              %(M0)f
        mu2             %(MU2)f
        perturbative    %(PERTURBATIVE)d
        n_burnin        %(N_BURNIN)d
        n_samples       %(N_SAMPLES)d
        ''' % {'M_LAT': M_lat,
               'T_FINAL':self.T_final,
               'M0':self.m0,
               'MU2':self.mu2,
               'PERTURBATIVE':self.perturbative,
               'N_BURNIN':self.n_burnin,
               'N_SAMPLES':self.n_samples}
        with open('parameters.in','w') as f:
            print >> f, parameters 
        output = subprocess.check_output([self.executable])
        for line in output.split('\n'):
            m = re.match('^ *mean *= *(.*) *variance *(.*)$',line)
            if m:
                mean = float(m.group(1))
                variance = float(m.group(2))
        return variance

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
        try:
            with open('variance.dat','r') as f:
                for line in f.readlines():
                    m = re.match(' *# *T_final *= *(.*)',line)
                    if m:
                        match_T_final = abs(self.T_final - float(m.group(1))) < tolerance
                    m = re.match(' *# *m0 *= *(.*)',line)
                    if m:
                        match_m0 = abs(self.m0 - float(m.group(1))) < tolerance
                    m = re.match(' *# *mu2 *= *(.*)',line)
                    if m:
                        match_mu2 = abs(self.mu2 - float(m.group(1))) < tolerance
        except:
            recreate_data = True
        recreate_data = not (match_T_final and match_m0 and match_mu2)
        if (recreate_data):
            with open('variance.dat','w') as f:
                print >> f, "# T_final = ",self.T_final
                print >> f, "# m0      = ",self.m0
                print >> f, "# mu2     = ",self.mu2
                for M_lat in self.M_lat_list:
                    variance = self._run_variance(M_lat)
                    print >> f, ('%12.8f' % M_lat)+'  '+('%12.8f' % variance)

    def plot_variance_decay(self):
        '''Plot variance of differences as a function of the lattice size
        '''
        self._run_all()
        M_lat = []
        variance = []
        with open('variance.dat') as f:
            for line in f.readlines():
                if (not (re.match(' *#.*',line))):
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
        var0 = 2.*self._singlelevel.Xsquared_continuum()**2
        p_var0 = plt.plot([M_lat[0],M_lat[-1]],[var0,var0],
                          linewidth=2,
                          linestyle='--')[0]
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
        self._C_deltaV = math.exp(a_0)
        p_fit = plt.plot(M_lat,self._C_deltaV*M_lat**a_1,
                         linewidth=2,
                         color='red')[0]

        plt.legend((p_var,p_var0,p_fit,p_lin),
                   (r'Var[$Y_\ell$]',
                    r'Var[$\langle X^2\rangle$]',
                    'linear fit',
                    r'$\propto M_\ell^{-1}$'),'center right')
        plt.savefig('variance_decay.pdf',bbox_inches='tight')

    @property
    def C_deltaV(self):
        if (self._C_deltaV == None):
            self.plot_variance_decay()
        return self._C_deltaV

class CostAnalysis(object):
    '''
    Class for comparative cost analysis

    :arg T_final: Final time
    :arg m0: Particle mass m_0
    :arg mu2: Quadratic potential parameter \mu^2
    '''
    def __init__(self,T_final,m0,mu2):
        self.T_final = T_final
        self.m0 = m0
        self.mu2 = mu2
        self.M0 = 8
        self.singlelevel = SingleLevel(self.T_final,self.m0,self.mu2)
        self.multilevel = MultiLevel(self.T_final,self.m0,self.mu2)

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

    print " m0      = ", m0
    print " mu^2    = ", mu2
    print " T_final = ", T_final
    print
    
    singlelevel = SingleLevel(T_final,m0,mu2)
    print " C_{bias} = ",singlelevel.Xsquared_bias_constant()
    print " V        = ",2.*singlelevel.Xsquared_continuum()**2
    singlelevel.plot_bias()

    multilevel = MultiLevel(T_final,m0,mu2)
    multilevel.plot_variance_decay()
    print "delta C_V = ",multilevel.C_deltaV

    cost_analysis = CostAnalysis(T_final,m0,mu2)
    cost_analysis.plot_cost()
    

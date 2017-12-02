'''Theoretical performance anAalysis of Monte Carlo methods'''
import math
import numpy as np
from sympy import *
from matplotlib import pyplot as plt

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
        return (self.Xsquared(M_ref)-self.Xsquared_continuum())*M_ref**2
    
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

    def plot_cost(self):
        '''Plot the cost as a function of the lattice spacing.

        If M is the number of time slices and N the number of MC samples, 
        then balancing the bias- and sampling-error gives a cost of

        Cost = N*M 
        '''
        M_lat = np.array([8,16,32,64,128,256])
        C_bias = abs(self.Xsquared_bias_constant())
        variance = 2.*self.Xsquared_continuum()
        epsilon = math.sqrt(2.)*C_bias*(1.0*M_lat)**(-2)
        cost = 2.**(5./4.)*variance*math.sqrt(C_bias)*epsilon**(-2.5)
        plt.clf()
        ax = plt.gca()
        ax.set_xscale('log')
        ax.set_yscale('log')
        p_bias = plt.plot(epsilon,cost,
                          linewidth=2,
                          color='blue',
                          marker='o',
                          markersize=8)[0]
        eps_ref = epsilon[-1]
        cost_ref = cost[-1]
        rho = 10.0
        p_theory = plt.plot([eps_ref,rho*eps_ref],
                            [0.25*cost_ref,0.25*rho**(-2.5)*cost_ref],
                            linewidth=2,
                            color='black')[0]
        ax.set_xlabel(r'tolerance $\epsilon$')
        ax.set_ylabel(r'Cost $N\cdot M$')
        plt.legend((p_bias,p_theory),
                   ('cost',r'$\propto \epsilon^{-2.5}$'),
                   'upper right')
        plt.savefig('cost_singlelevel.pdf',bbox_inches='tight')


        pass
        
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
    singlelevel.plot_cost()

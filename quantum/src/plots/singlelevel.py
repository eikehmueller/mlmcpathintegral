'''Cost parameters of single level method'''
import math
import numpy as np
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

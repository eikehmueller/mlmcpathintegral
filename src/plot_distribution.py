import re
import sys
import math
import numpy as np
from matplotlib import pyplot as plt

'''
  Script for visualising data generated by the test_distribution.cc
  main program.

  This script processes text files in the format given below. It plots
  a binned distribution and compares the evaluated distribution to the
  exact analytical expression. The file header depends on the considered 
  distribution and is parsed to extract the relevant parameters.

  -------------------------------------------------
  BesselProductDistribution
    beta = 4
    x_p  = 2.82743
    x_m  = 0
 
  n_samples = 1000
  n_points = 129
 
  ==== samples ====
  2.05174
  1.82559
  0.806089
  [...]
  1.50355
 
  ==== points ====
  -3.14159 0.0371872
  -3.09251 0.0364034
  [...]
  3.09251 0.0385538
  3.14159 0.0371872
  -------------------------------------------------
'''

class Distribution(object):
    '''Class for wrapping a distribution

    :arg samples: Samples drawn for distribution
    :arg X: x-coordinates of points to be plotted
    :arg Y: y-coordinates of points to be plotted
    '''
    def __init__(self,samples,X,Y):
        self.samples = samples
        self.X = X
        self.Y = Y
        self.label = 'UNKNOWN'

    '''Plot a histogram of the samples and the distribution 
       and save to file

    :arg filename: Name of file to write to
    :arg bins: Number of bins in sample histogram
    '''
    def plot(self,filename,bins=64):
        plt.clf()
        X_exact = np.arange(-np.pi,np.pi,1.E-3)
        Y_exact = np.vectorize(self.f_exact)(X_exact)
        plt.hist(self.samples,
                     histtype='stepfilled',
                     alpha=0.5,
                     density=True,
                     bins=bins,
                     label='samples')
        plt.plot(self.X,self.Y,
                     linewidth=2,
                     linestyle='-',
                     color='black',
                     label='distribution')
        plt.plot(X_exact,Y_exact,
                     linewidth=2,
                     linestyle='--',
                     color='red',
                     label='exact')
        ax = plt.gca()
        ax.set_xlim(-1.05*np.pi,1.05*np.pi)
        ax.set_xticks((-np.pi,-0.5*np.pi,0,0.5*np.pi,np.pi))
        ax.set_xticklabels((r'$-\pi$',r'$-\frac{\pi}{2}$','$0$',r'$\frac{\pi}{2}$',r'$\pi$'))
        plt.legend(loc='upper right')
        plt.title(self.label)
        plt.savefig(filename,bbox_inches='tight')


class ExpSin2Distribution(Distribution):
    '''Wrapper for ExpSin2Distribution

    :arg samples: Samples drawn for distribution
    :arg X: x-coordinates of points to be plotted
    :arg Y: y-coordinates of points to be plotted
    :arg beta: parameter beta
    '''
    def __init__(self,samples,X,Y,sigma):
        super().__init__(samples,X,Y)
        self.label = 'ExpSin2Distribution'
        self.sigma = sigma
        self.Znorm_inv = np.exp(0.5*self.sigma)/(2.*np.pi*np.i0(0.5*self.sigma))

    '''Exact value of distribution at a given point
    
    :arg x: Point at which the distribution is evaluated
    '''
    def f_exact(self,x):
        return self.Znorm_inv*np.exp(-self.sigma*np.sin(0.5*x)**2)

class ExpCosDistribution(Distribution):
    def __init__(self,samples,X,Y,beta,x_p,x_m):
        super().__init__(samples,X,Y)
        self.label = 'ExpCosDistribution'
        self.beta = beta
        self.x_p = x_p
        self.x_m = x_m
        self.Znorm = 2.*np.pi*np.i0(2.*beta*np.cos(0.5*(self.x_p-self.x_m)))

    '''Exact value of distribution at a given point
    
    :arg x: Point at which the distribution is evaluated
    '''
    def f_exact(self,x):
        return 1./self.Znorm*np.exp(self.beta*(np.cos(x-self.x_p)+np.cos(x-self.x_m)))
    
class BesselProductDistribution(Distribution):
    def __init__(self,samples,X,Y,beta,x_p,x_m):
        super().__init__(samples,X,Y)
        self.label = 'BesselProductDistribution'
        self.beta = beta
        self.x_p = x_p
        self.x_m = x_m
        # Number of expansion terms used for evaluation of
        # normalisation constant
        self.kmax = 16
        # Maximal number of terms used in the double sums that
        # define the terms A_k in the normalisation constant
        # expansion
        self.n_max = 32
        # Compute normalisation constant
        self.Znorm = 0
        # A_k
        for k in range(0,self.kmax+1):
            s = 0.0
            for n in range(k,self.n_max):
                for m in range(k,self.n_max):
                    s += (0.5*self.beta)**(2*(n+m))*math.comb(2*n,n-k)*math.comb(2*m,m-k)/(math.factorial(n)*math.factorial(m))**2
            if (k==0):
                A_k = 2*math.pi*s
            else:
                A_k = 4*math.pi*s
            self.Znorm += A_k*math.cos(k*(self.x_p-self.x_m))

    '''Exact value of distribution at a given point
    
    :arg x: Point at which the distribution is evaluated
    '''
    def f_exact(self,x):
        return 1./self.Znorm*np.i0(2*self.beta*np.cos(0.5*(x-self.x_p)))*np.i0(2*self.beta*np.cos(0.5*(x-self.x_m)))

'''Read data in the format given above from a text file
   Returns a distribution wrapper of the type saved in the file

:arg: filename Name of file to read
'''
def read_data(filename):
    mode = None
    j = 0 # Counter for points
    k = 0 # Counter for samples
    param = {} # Dictionary for saving distribution parameters
    with open(filename) as f:
        for line in f.readlines():
            # Read in samples
            if (mode == 'ParseSamples'):
                if (j<n_samples) and (line.strip() != ''):
                    samples[j] = float(line)
                    j += 1
            elif (mode == 'ParsePoints'):
                if (k<n_points) and (line.strip() != ''):
                    x,y = line.split()
                    X[k] = float(x)
                    Y[k] = float(y)
                    k += 1
            else:
                # Work out which name of distribution
                m = re.match(' *(.*Distribution)',line)
                if m:
                    distribution = m.group(1)
                # Work out number of samples contained in file
                m = re.match(' *n_samples *= *([0-9]+)',line)
                if m:
                    n_samples = int(m.group(1))
                # Work out number of points contained in file
                m = re.match(' *n_points *= *([0-9]+)',line)
                if m:
                    n_points = int(m.group(1))
                # Parse any other parameters
                m = re.match(' *([a-zA-Z0-9_]+) *= *([0-9\.\-\+Ee]+)',line)
                if m:
                    param[m.group(1)] = m.group(2)
            # Check for line marking start of samples
            m = re.match('.*==== samples ====.*',line)
            if m:
                mode = 'ParseSamples'
                samples = np.zeros(n_samples)
            m = re.match('.*==== points ====.*',line)
            if m:
                mode = 'ParsePoints'
                X = np.zeros(n_points)
                Y = np.zeros(n_points)
    print ('distribution = ',distribution)
    print ('n_samples = ',n_samples)
    print ('n_points = ',n_points)
    print ('distribution parameters:')
    for (key,value) in param.items():
        if not ((key == 'n_samples') or (key == 'n_points')):
            print ('  ',key,' = ',value)
    if (distribution == 'ExpSin2Distribution'):
        return ExpSin2Distribution(samples,X,Y,float(param['sigma']))
    elif (distribution == 'ExpCosDistribution'):
        return ExpCosDistribution(samples,X,Y,
                                      float(param['beta']),
                                      float(param['x_p']),
                                      float(param['x_m']))
    elif (distribution == 'BesselProductDistribution'):
        return BesselProductDistribution(samples,X,Y,
                                             float(param['beta']),
                                             float(param['x_p']),
                                             float(param['x_m']))
    else:
        print ('ERROR: Unknown distribution: ',distribution)
        
#################################################################
###                         M A I N                           ###
#################################################################
if (__name__ == '__main__'):
    nbins = 64 # Number of bins used for sample histogram
    if (len(sys.argv) != 2):
        print ('Usage: python3 '+sys.argv[0]+' FILENAME')
        sys.exit(-1)
    filename = sys.argv[1]
    distribution = read_data(filename)
    distribution.plot('distribution.pdf')
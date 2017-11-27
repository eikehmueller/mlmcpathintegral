import numpy as np
from sampler import *

class Action(object):
    '''Base class for action of a one dimensional quantum mechanical
    particles described by its path vector (x_0,x_1,\dots,x_{M-1})
    where x_j = x(a*j) and a = T/M. Periodic boundary conditions x_M = x_0 
    are assumed.

    :arg M: Number of lattice points = timesteps
    :arg T: Final timestep
    :arg m0: Mass m_0
    '''
    def __init__(self,M,T,m0=1.0):
        self.M = M
        self.T = T
        self.a = T/M
        self.m0 = m0

    def evaluate(self,x):
        '''Evaluate action S[x] on a path x = (x_0,x_1,...,x_{M-1})

        :arg x: Path vector, has to be of length M
        '''
        assert(len(x)==self.M)

class HarmonicOscillatorAction(Action,Sampler):
    '''Action for harmonic oscillator with mass m_0 and 
    potential 1/2*m_0*\mu^2*x^2.

    Note that since we can sample directly in this case, we also derive
    from the Sampler base class.

    :arg M: Number of lattice points = timesteps
    :arg T: Final timestep
    :arg m0: Mass m_0
    :arg mu2: Parameter \mu^2
    '''

    def __init__(self,M,T,m0=1.0,mu2=1.0):
        super(HarmonicOscillatorAction,self).__init__(M,T,m0)
        self.mu2 = mu2
        # Mean of Gaussian distribution
        self._mean = np.zeros(self.M)
        self._Sigma = self._covariance()

    def evaluate(self,x):
        '''Evaluate action S[x] on a path x = (x_0,x_1,...,x_{M-1})

        :arg x: Path vector, has to be of length M
        '''        
        assert(len(x)==self.M)
        S = 0.0
        ainv2 = 1./self.a**2
        for j in range(self.M):
            S += (ainv2*(x[j]-x[j-1])**2+self.mu2*x[j]**2)
        S *= 0.5*self.a*self.m0
        return S
       
    def draw(self):
        '''Generate a sample x from the Gaussian distribution induced
        by the action, i.e. from C*exp(-1/2*x^T.\Sigma^{-1}.x)
        where the matrix \Sigma is defined by \Sigma = G^{-1}
        and -1/2*x^T.G.x = S[x]'''
        return np.random.multivariate_normal(self._mean,self._Sigma)

    def x2continuum(self):
        '''Return <x^2> in the continuum (but not infinite volume) limit'''
        E = math.exp(-self.T*math.sqrt(self.mu2))
        return 0.5/math.sqrt(self.mu2)*(1.+E)/(1.-E)

    def _covariance(self):
        '''Construct and return the covariance
        matrix of the multivariate Gaussian'''
        # Construct matrix G
        G = np.zeros((self.M,self.M))
        d_tmp = self.a*self.m0*self.mu2 + 2.*self.m0/self.a
        c_tmp = -self.m0/self.a
        for j in range(self.M):
            G[j,j] = d_tmp
            G[j,(j+1)%self.M] = c_tmp
            G[j,(j-1)%self.M] = c_tmp
        self._Sigma = np.linalg.inv(G)
        return self._Sigma

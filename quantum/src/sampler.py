import math
import numpy as np

class Sampler(object):
    '''Base class for sampler'''
    def __init__(self):
        pass

    def sample(self):
        pass

class TwoLevelMetropolisSampler(object):
    '''Two-level sampler. Given a way of generating coarse
    samples theta_{\ell-1} from a coarse_sampler, create
    fine samples theta_{\ell} using a two-level Metropolis-Hastings
    process.

    Draw fine state from free distribution conditioned by
    start-point x_- and end-point x_+
    with
    
    p(x) ~ exp[-m_0/a (x-0.5*(x_+ + x_-))^2]
        
    :arg coarse_sampler: Sampler for coarse distribution
    :arg coarse_action: Action on coarse level
    :arg fine_action: Action on fine level
    :arg verbosity: Verbosity level (0=quiet, 1=print)
    '''
    def __init__(self,
                 coarse_sampler,
                 coarse_action,
                 fine_action,
                 verbosity=0):
        self._coarse_sampler = coarse_sampler
        self._coarse_action = coarse_action
        self._fine_action = fine_action
        self._verbosity = verbosity
        self.M = self._fine_action.M
        self.delta = math.sqrt(self._fine_action.a)
        # Initial fine state
        self.theta_fine = np.zeros(self.M)
        self._free_width = math.sqrt(0.5*self._fine_action.a/self._fine_action.m0)

    def _conditioned_free_action(self,theta):
        '''Calculate action 
        
        S = m_0/a*\sum_{j=0}^{M/2}(theta_{2j}-0.5*(theta_{2j}+theta_{2j+1}))^2
        
        which gives the conditioned distribution of the free particles at the 
        fine sites,

        q_{ML}^{\ell,F}(\theta_{\ell,F}|\theta_{\ell,C})

        where \theta_{\ell} = [\theta_{\ell,C},\theta_{\ell,F}]
       
        :arg theta: State vector
        '''
        S = 0.0
        for j in range(self.M/2):
            S += (theta[2*j+1] - 0.5*(theta[2*j]+theta[(2*j+2)%self.M]))**2
        return self._fine_action.m0/self._fine_action.a*S

    def draw(self):
        '''Create new correlated fine- and coarse-level state pair.

        Draw a sample from the coarse level distribution.
        Then construct new fine level trial state
        \theta'_{\ell} = \theta'_{C,\ell}, \theta'_{F,\ell}
        By setting taking the coarse modes from the next-coarser
        level (\theta'_{\ell,C} = \Theta_{\ell-1}^{n+1} and filling 
        in the fine model by drawing from the free action, conditioned
        on the coarse modes.

        Finally, accept/reject to ensure that the fine level state has the
        correct distribution.
        
        List of variables                          length
           theta_coarse    = \Theta_{\ell}^{n+1}   M/2
           theta_prime     = \theta'_{\ell}        M
           self.theta_fine = \theta_{\ell}^{n}     M
        '''

        # Obtain new coarse level sample \Theta_{\ell-1}^{n+1}
        theta_coarse = self._coarse_sampler.sample()
        # Populate fine level trial state \theta'_{\ell}
        theta_prime = self.theta_fine.copy()
        for j in range(self.M/2):
            # coarse modes (C)
            theta_prime[2*j] = theta_coarse[j]
            # fine modes (F)
            mu = 0.5*(theta_coarse[j]+theta_coarse[(j+1)%(self.M/2)])
            theta_prime[2*j+1] = np.random.normal(mu,self._free_width)
        # Calculate the difference in level-\ell actions,
        # required to calculate the ratio
        #   \pi^{\ell}(\theta'_\ell)/\pi^{\ell}(\theta_\ell^{n})
        deltaS_fine = self._fine_action.evaluate(theta_prime) \
                    - self._fine_action.evaluate(self.theta_fine)
        # Calculate the difference in level-(\ell-1) actions,
        # required to calculate the ratio
        #   \pi^{\ell-1}(\theta_{\ell,C}^{n})/\pi^{\ell-1}(\theta'_{\ell,C})
        deltaS_coarse = self._coarse_action.evaluate(self.theta_fine[::2]) \
                      - self._coarse_action.evaluate(theta_coarse)
        # Calculate the difference in free level-\ell actions,
        # required to calculate the ratio
        #   q_{ML}^{\ell,F}(\theta_{\ell,F}^{n}|\theta_{\ell,C}^{n}) / 
        #   q_{ML}^{\ell,F}(\theta'_{\ell,F}|\theta'_{\ell,C})
        deltaS_trial = self._conditioned_free_action(self.theta_fine) \
                     - self._conditioned_free_action(theta_prime)
        deltaS = deltaS_fine + deltaS_coarse + deltaS_trial
        if (self._verbosity > 0):
            print 'delta S [fine]   = ', deltaS_fine
            print 'delta S [coarse] = ', deltaS_coarse
            print 'delta S [trial]  = ', deltaS_trial
            print 'delta S          = ', deltaS
            print 'accepted         = ', accept
            print
        accept = False
        if (deltaS < 0.0):
            accept = True
        else:
            t = math.exp(-deltaS)
            accept = (np.random.random()<t)
        if (accept):        
            self.theta_fine = theta_prime
        return (theta_coarse,self.theta_fine)

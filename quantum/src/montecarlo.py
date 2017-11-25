import math
import numpy as np
from sampler import *

class MonteCarloSingleLevel(object):
    '''Single level MC Monte Carlo sampler
    
    :arg action: action to use
    :arg sample: Sampler object
    :arg qoi: Quantity of interest
    :arg n_samples: Number of samples for evaluation of QoI
    :arg n_burnin: Number of samples for burnin
    '''
    def __init__(self,action,
                 sampler,
                 qoi,
                 n_samples,
                 n_burnin):
        self._action = action
        self._sampler = sampler
        self._qoi = qoi
        self._n_samples = n_samples
        self._n_burnin = n_burnin

    def evaluate(self):
        '''Calculate the MCMC estimate (and its error) of the QoI
        by burning in and then evaluating the QoI on n_samples samples
        '''
        for i in range(self._n_burnin):
            self._sampler.draw()
        S = 0.0
        S_sq = 0.0
        for i in range(self._n_samples):
            x = self._sampler.draw()
            qoi_tmp = self._qoi(x)
            S += qoi_tmp
            S_sq += qoi_tmp**2

        mean = S/self._n_samples
        variance = 1./(self._n_samples*(self._n_samples-1))*(S_sq - self._n_samples*mean**2)
        return mean, math.sqrt(variance)

class MonteCarloTwoLevel(object):
    '''Two level MC Monte Carlo sampler
    
    :arg coarse_sampler: Exact sampler for coarse level
    :arg coarse_action: Coarse level action
    :arg action: action to use
    :arg sample: Sampler object
    :arg qoi: Quantity of interest
    :arg n_samples: Number of samples for evaluation of QoI
    :arg n_burnin: Number of samples for burnin
    '''
    def __init__(self,
                 coarse_sampler,
                 coarse_action,
                 action,
                 qoi,
                 n_samples,
                 n_burnin):
        self._qoi = qoi
        self._n_samples = n_samples
        self._n_burnin = n_burnin
        self._twolevel_sampler = TwoLevelMetropolisSampler(coarse_sampler,
                                                           coarse_action,
                                                           action,
                                                           verbosity=0)

    def evaluate(self):
        '''Calculate the MCMC estimate (and its error) of the QoI
        by burning in and then evaluating the QoI on n_samples samples
        '''
        for i in range(self._n_burnin):
            x_coarse, x_fine = self._twolevel_sampler.draw()
        S = 0.0
        S_sq = 0.0
        for i in range(self._n_samples):
            x_coarse, x_fine = self._twolevel_sampler.draw()
            qoi_tmp = self._qoi(x_fine)
            S += qoi_tmp
            S_sq += qoi_tmp**2

        mean = S/self._n_samples
        variance = 1./(self._n_samples*(self._n_samples-1))*(S_sq - self._n_samples*mean**2)
        return mean, math.sqrt(variance)

    def evaluate_difference(self):
        '''Calculate the difference between the QoI on the fine- and coarse
        level and return this difference and its variance.
        '''
        for i in range(self._n_burnin):
            x_coarse, x_fine = self._twolevel_sampler.draw()
        S = 0.0
        S_sq = 0.0
        for i in range(self._n_samples):
            x_coarse, x_fine = self._twolevel_sampler.draw()
            qoi_tmp = self._qoi(x_fine)-self._qoi(x_coarse)
            S += qoi_tmp
            S_sq += qoi_tmp**2

        mean = S/self._n_samples
        variance = 1./(self._n_samples-1)*(S_sq - self._n_samples*mean**2)
        return mean, variance

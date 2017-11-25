import numpy as np
import math
from action import *
from montecarlo import *
from quantityofinterest import *

M = 64
T = 1.0
m0 = 1.0
mu2 = 1.0
n_burnin = 10000
n_samples = 1000000

qoi = QoIXsquared()
action = HarmonicOscillatorAction(M,T,m0,mu2)
coarse_action = HarmonicOscillatorAction(M/2,T,m0,mu2)

# Print out parameters

print 'Number of timesteps = ', M
print 'Final time T        = ', T
print 'mass m_0            = ', m0
print 'mu^2                = ', mu2
print '# burn-in samples   = ', n_burnin
print '# samples           = ', n_samples
print

print 'Continuum limit <x^2> = ',action.x2continuum()
print

montecarlo_singlelevel = MonteCarloSingleLevel(action,
                                               action,
                                               qoi,
                                               n_samples,
                                               n_burnin)

x_sq_mean, x_sq_error = montecarlo_singlelevel.evaluate()

print '*** single-level method (direct sampling) ***'
print '<x^2> = ' + ('%8.4f' % x_sq_mean) + ' +/- ' + ('%8.4f' % x_sq_error)
print

montecarlo_twolevel = MonteCarloTwoLevel(coarse_action,
                                         coarse_action,
                                         action,
                                         qoi,
                                         n_samples,
                                         n_burnin)
x_sq_mean, x_sq_error = montecarlo_twolevel.evaluate()

print '*** two-level method ***'
print '<x^2> = ' + ('%8.4f' % x_sq_mean) + ' +/- ' + ('%8.4f' % x_sq_error)
print

# Variance decay:
with open('variance.dat','w') as f:
    print >> f, '# a   Var(Y)'
    for logM in range(3,8):
        M = 2**logM
        action = HarmonicOscillatorAction(M,T,m0,mu2)
        coarse_action = HarmonicOscillatorAction(M/2,T,m0,mu2)
        montecarlo_twolevel = MonteCarloTwoLevel(coarse_action,
                                                 coarse_action,
                                                 action,
                                                 qoi,
                                                 n_samples,
                                                 n_burnin)
        diff_mean, diff_var = montecarlo_twolevel.evaluate_difference()
        print ('%4d' % M)+' '+('%8.4f' % diff_var)
        print >> f,('%12.8e' % (T/M))+' '+('%12.8f' % diff_var)

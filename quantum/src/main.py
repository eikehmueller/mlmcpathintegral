import numpy as np
import math
from action import *

N = 32
T = 1.0
m0 = 1.0
mu2 = 1.0
n_burnin = 1000
n_samples = 10000

action = HarmonicOscillatorAction(N,T,m0,mu2)
coarse_action = HarmonicOscillatorAction(N/2,T,m0,mu2)

# Print out parameters

print 'Number of timesteps = ', N
print 'Final time T        = ', T
print 'mass m_0            = ', m0
print 'mu^2                = ', mu2
print '# burn-in samples   = ', n_burnin
print '# samples           = ', n_samples
print

print 'Continuum limit <x^2> = ',action.x2continuum()
print

# Sample with single-level method
x_sq = 0.0
for i in range(n_samples):
    x = action.sample()
    x_sq += 1./N*np.dot(x,x)

x_sq /= n_samples

print '*** single-level method (direct sampling) ***'
print '<x^2> = ', x_sq
print

# Sample with two-level method
twolevel_sampler = TwoLevelMetropolisSampler(coarse_action,
                                             coarse_action,
                                             action,
                                             verbosity=0)
x_sq = 0.0
for i in range(n_burnin):
    x_coarse, x_fine = twolevel_sampler.generate()
for i in range(n_samples):
    x_coarse, x_fine = twolevel_sampler.generate()
    x_sq += 1./N*np.dot(x_fine,x_fine)

x_sq /= n_samples

print '*** two-level method ***'
print '<x^2> = ', x_sq
print

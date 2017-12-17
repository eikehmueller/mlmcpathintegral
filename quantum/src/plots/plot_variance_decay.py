'''Theoretical performance Analysis of Monte Carlo methods'''
import sys
import subprocess
import math
import numpy as np
from multilevel import *

    
if (__name__ == '__main__'):
    # Numerical values for parameters
    T_final = 1.0
    m0 = 1.0
    mu2 = -1.0
    lmbda = 1.0

    action = 'quartic oscillator'

    print " action  = ", action
    print " m0      = ", m0
    print " mu^2    = ", mu2
    print " lambda  = ", lmbda
    print " T_final = ", T_final
    print
    
    multilevel = MultiLevel(action,T_final,m0,mu2,lmbda)
    multilevel.plot_variance_decay()
    print "delta C_V = ",multilevel.C_deltaV

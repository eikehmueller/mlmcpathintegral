###########################################
# General parameters
###########################################
general:
  do_singlelevelmc = false # Run single level MC algorithm?
  do_twolevelmc = true # Run two level MC algorithm?
  do_multilevelmc = false # Run multilevel MC algorithm?  
  action = 'rotor' # action to use. Allowed values are:
                                  # [harmonicoscillator,
                                  #  quarticoscillator,
                                  #               doublewell,
                                  #  rotor]

###########################################
# Lattice parameters
###########################################
lattice:
  M_lat = %(MLAT)d        # Number of lattice sites
  T_final = %(TFINAL)f     # Final time

###########################################
# Statistics parameters
###########################################
statistics:
  n_autocorr_window = 20   # Size of window over which to measure
                           # autocorrelations
  n_min_samples_corr = 1000 # Minimal number of samples for correlated
                           # estimators
  n_min_samples_qoi = 1000  # Minimal number of samples for uncorrelated
                           # estimators

###########################################
# Harmonic oscillator action parameters
###########################################
harmonicoscillator:
  m0 = 1.0           # Unrenormalised mass
  mu2 = 1.0          # Curvature of potential
  renormalisation = 'none' # Renormalisation [none, perturbative, exact]

###########################################
# Double well action parameters
###########################################
doublewell:
  m0 = 1.0           # Unrenormalised mass
  mu2 = 1.0          # Curvature of potential
  lambda = 1.0       # Parameter lambda
  sigma = 1.0        # Parameter sigma

###########################################
# Quartic oscillator action parameters
###########################################
quarticoscillator:
  m0 = 1.0           # Unrenormalised mass
  mu2 = -1.0          # Curvature of potential
  lambda = 1.0       # Parameter lambda

###########################################
# Rotor action parameters
###########################################
rotor:
  m0 = %(M0)f          # Unrenormalised mass

###########################################
# Single level Monte Carlo parameters
###########################################
singlelevelmc:
  n_burnin = 1000     # Number of burnin samples
  epsilon = 1.E-3
  sampler = 'cluster'    # Sampler to use [HMC, cluster, exact]

###########################################
# Two level Monte Carlo parameters
###########################################
twolevelmc:
  n_burnin = 1000      # Number of burnin samples
  n_samples = %(NSAMPLES)d    # Number of samples
  coarsesampler = 'cluster' # Sampler to use [HMC, cluster, exact]

###########################################
# Multilevel Monte Carlo parameters
###########################################
multilevelmc:
  n_level = 3      # Number of levels
  n_burnin = 100      # Number of burnin samples
  epsilon = 0.000137
  coarsesampler = 'exact' # Sampler to use [HMC, cluster, exact]
  show_detailed_stats = true

###########################################
# HMC sampler parameters
###########################################
hmc:
  T = 1.0        # HMC integration time
  dt = 0.01       # HMC time step
  n_burnin = 100  # Number of burnin samples

###########################################
# Cluster sampler parameters
###########################################
clusteralgorithm:
  n_burnin = 1000  # Number of burnin samples
  n_updates = 10  # Number of cluster updates between steps

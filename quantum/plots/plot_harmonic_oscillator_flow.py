import numpy as np
import math
from matplotlib import pyplot as plt

import matplotlib as mpl

mpl.rcParams.update({'font.size': 22})

def coarsen(a,m0,mu_sq):
    m0_c = m0/(1.+a**2*mu_sq/(2.*m0))
    mu_sq_c = mu_sq*(1.+a**2*mu_sq/(4.*m0))/(1.+a**2*mu_sq/(2.*m0))
    return 2*a, m0_c, mu_sq_c

n = 8

a_min = 2.**(1-n)
a = a_min 
m0_star = 0.01
m0 = m0_star 
mu_sq_star = 1.0
mu_sq = mu_sq_star

a_list = []
m0_list = []
mu_sq_list = []

for i in range(n):
    a_list.append(a)
    m0_list.append(m0)
    mu_sq_list.append(mu_sq)
    a_max = a
    a, m0, mu_sq = coarsen(a,m0,mu_sq)

a = np.array(a_list)
m0 = np.array(m0_list)
mu_sq = np.array(mu_sq_list)

plt.clf()
ax = plt.gca()
ax.set_yscale('log')
ax.set_ylim(1.E-6,10.)
ax.set_xlim(1./2**n,2.)
ax.set_xscale('log')
p_m0 = plt.plot(a,m0,linewidth=2,color='blue',marker='o',markersize=10,markeredgecolor='blue')[0]
plt.plot([a_min,a_max],[m0_star,m0_star],linewidth=2,color='blue',linestyle='--')
ax.annotate('$m_0(a^*)$', xy=(math.sqrt(a_min*a_max), 2.*m0_star),color='blue')
p_mu_sq = plt.plot(a,mu_sq,linewidth=2,color='red',marker='s',markersize=10,markeredgecolor='red')[0]
plt.plot([a_min,a_max],[mu_sq_star,mu_sq_star],linewidth=2,color='red',linestyle='--')
ax.annotate('$\mu^2(a^*)$', xy=(math.sqrt(a_min*a_max), 2.*mu_sq_star),color='red')

plt.legend((p_m0,p_mu_sq),('$m_0(a)$','$\mu^2(a)$'),'lower left')
ax.set_xlabel('a')
ax.set_ylabel('$m_0$ or $\mu^2$')
plt.savefig('harmonic_oscillator_flow.pdf',bbox_inches='tight')

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
mpl.rcParams.update({'font.size': 22})

a_list = []
varY_list = []
with open('variance.dat') as f:
    for line in f.readlines():
        try:
            a, varY = line.split()
            a_list.append(float(a)), varY_list.append(float(varY))
        except:
            pass

a_list = np.array(a_list)
varY_list = np.array(varY_list)
plt.clf()
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'lattice spacing $a=2^{-\ell}T$')
ax.set_ylabel(r'$\operatorname{Var}(Y^\ell)$')
p_var = plt.plot(a_list,varY_list,
                 linewidth=2,color='blue',marker='o',
                 markersize=8)[0]

p_linear = plt.plot([0.02,0.16],[0.02,0.16],linewidth=2,color='black')[0]

plt.legend((p_var,p_linear),('variance','linear decay'),'upper left')

plt.savefig('variance_decay.pdf',bbox_inches='tight')

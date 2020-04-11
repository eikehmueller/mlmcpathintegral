import numpy as np
import pandas
from matplotlib import pyplot as plt

def c(data):
    avg = np.average(data)
    n_array = data.size-np.arange(data.size)
    tmp = np.correlate(data,data,mode='full')
    result = tmp[tmp.size//2:]*1./n_array-avg**2
    return result/result[0]

def tau_int(data_c,M):
    n = data_c.size
    tmp = 0.0
    for k in range(1,M):
        tmp += (1.-k/n)*data_c[k]
    return 1.+2*tmp

# Create a dtype with the binary data format and the desired column names
dt = np.dtype([('qoi', 'f8')])
filename='qoi_M128.dat'
data = np.fromfile(filename, dtype=dt)
df = pandas.DataFrame.from_records(data)

X = np.asarray(df.index)
Y = np.asarray(df['qoi'])
acY = c(Y)
plt.clf()
nmax=2000
tauint = tau_int(acY,1000)
print ('tau_{int} = ',tauint)
plt.plot(acY[:nmax],linewidth=2,color='blue',marker='o',markersize=4,markeredgewidth=2,markerfacecolor='white')
X = np.arange(0,nmax,1.E-2)
plt.plot(X,np.exp(-X/(tauint/2)),linewidth=2,color='red')
plt.plot([0,nmax],[0,0],linewidth=2,color='black',linestyle='--')
plt.savefig('qoi.pdf',bbox_inches='tight')

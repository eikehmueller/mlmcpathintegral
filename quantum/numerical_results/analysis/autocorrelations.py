import numpy as np
import pandas
from matplotlib import pyplot as plt
import math
                 
class AutoCorrelation(object):

    def __init__(self,data):
        '''Class for analysing the Autocorrelations of a time series
        
        :arg data: time-series data
        '''
        self._data = data
        self._data_c = None
    
    def _c(self):
        '''Compute normalised autocorrelation function in temporal domain

        Let the (empirical) covariance function be

        \gamma(k) = 1/(N-k)\sum_{j=0}^{N-k-1} (X(j) - \mu)*(X(j+k) - \mu)

        where \mu = 1/N\sum_{j=0}^{N-1} X(j) is the mean.
        
        Then this method returns [c(0),c(1),...,c(N-1)] with 
        c(k) = \gamma(k)/\gamma(0)
        '''
        if (self._data_c) is None:
            avg = np.average(self._data)
            n_array = self._data.size-np.arange(self._data.size)
            tmp = np.correlate(self._data,self._data,mode='full')
            result = tmp[tmp.size//2:]*1./n_array-avg**2
            self._data_c = result/result[0]
        return self._data_c

    def plot(self,filename):
        data_c = self._c()
        tauint = self.tau_int()
        M = math.ceil(5*tauint)
        tauint_naive = self.tau_int_naive(M)
        print ('tau_{int}         = '+('%10.2f' % tauint))
        print ('tau_{int} [naive] = '+('%10.2f' % tauint_naive))
        plt.clf()
        ax = plt.gca()
        ax.set_xlabel('Step $k$')
        ax.set_ylabel('$c(k)$')
        plt.plot(data_c[:M],
                 linewidth=2,
                 color='blue',
                 marker='o',
                 markersize=4,
                 markeredgewidth=2,
                 markerfacecolor='white',
                 label='data')
        X = np.arange(0,M,1.E-2)
        plt.plot(X,np.exp(-X/(0.5*(tauint+1))),
                 linewidth=2,
                 color='red',
                 label=r'$exp^{-j/\tau_{\textrm{exp}}}$')
        plt.plot(X,np.exp(-X/(0.5*(tauint_naive+1))),
                 linewidth=2,
                 linestyle='--',
                 color='red',
                 label=r'$exp^{-j/\tau_{\textrm{exp}}^{(\textrm{naive})}}$')
        plt.plot([0,M],[0,0],linewidth=2,color='black',linestyle='--')
        plt.savefig(filename,bbox_inches='tight')
        
    def _C1(self,K,degree):
        ''' The function C_1(K,d) defined in the appendix of
        Heidelberger & Welch

        :arg K: Number of points in spectral space to use
        :arg degree: Polynomial degree for fit
        '''
        N = len(self._data)
        f_n = (4.*np.arange(1,K+1)-1.)/(2.*N)
        X = np.zeros((K,degree+1))
        for k in range(K):
            for j in range(degree+1):
                X[k,j] = f_n[k]**j
        sigma2 = 0.645*(np.linalg.inv(np.transpose(X)@X)[0,0])
        return np.exp(-0.5*sigma2)
    
    def tau_int_naive(self,M):
        '''
        Compute naive integrated autocorrelation time in temporal domain
        
        This uses the 'naive' definition of the integrated autocorrelation
        time, summing up to an upper value of M.

        \tau_{int} = 1 + 2*\sum_{k=1}^{M} (1 - k/N)*c(k)

        where c(k) is the function returned by the method _c
        
        :arg M: Size of window
        
        '''
        data_c = self._c()
        n = data_c.size
        tmp = 0.0
        for k in range(1,M):
            tmp += (1.-k/n)*data_c[k]
        return 1.+2.*tmp

    def tau_int(self,K=50,degree=2):
        return self._tau_int(self._data,K,degree)

    def _tau_int(self,data,K=50,degree=2):
        '''Computes the integrated autocorrelation time using the
        method in the paper

        Heidelberger, P. and Welch, P.D., 1981. A spectral method for
        confidence interval generation and run length control in simulations.
        Communications of the ACM, 24(4), pp.233-245.

        :arg K: Number of points in spectral space to use
        :arg degree: Degree of fit polynomial
        '''
        N = len(data)
        I_freq = np.abs(np.fft.fft(data)[:N//2])**2/N
        J = np.zeros(N//4)
        for n in range(N//4):
            J[n] = np.log(0.5*(I_freq[2*n-1]+I_freq[2*n]))
        X_fit = (4.*np.arange(1,K+1)-1.)/(2.*N)
        Y_fit = J[1:K+1]+0.270

        a_fit=np.polyfit(X_fit,Y_fit,deg=degree)
        var = np.var(data)
        p0 = self._C1(K,degree)*np.exp(a_fit[-1])
        return p0/var

def generate_celerite_data(n_data):

    np.random.seed(123457)

    # Build the celerite model:
    import celerite
    from celerite import terms
    kernel = terms.RealTerm(log_a=0.0, log_c=-6.0)
    kernel += terms.RealTerm(log_a=0.0, log_c=-2.0)

    # The true autocorrelation time can be calculated analytically:
    true_tau = sum(2*np.exp(t.log_a-t.log_c) for t in kernel.terms)
    true_tau /= sum(np.exp(t.log_a) for t in kernel.terms)
    true_tau_int = 2*true_tau - 1
    print ('True tau_{int} = ',true_tau_int)
    # Simulate a set of chains:
    gp = celerite.GP(kernel)
    t = np.arange(n_data)
    gp.compute(t)
    y = gp.sample(size=32)
    
    # Let's plot a little segment with a few samples:
    plt.plot(y[:3, :300].T)
    plt.xlim(0, 300)
    plt.xlabel("step number")
    plt.ylabel("$f$")
    plt.title(r"$\tau_{\mathrm{int}}^{(\mathrm{true})} = "+"{0:.0f}$".format(true_tau_int), fontsize=14);
    plt.savefig('celerite_process.pdf',bbox_inches='tight')
    print ('saved to file celerite_process.pdf')
    return y[0].T

def read_data(filename):
    print ('Reading data from file '+filename)
    # Create a dtype with the binary data format and the desired column names
    dt = np.dtype([('qoi', 'f8')])
    data = np.fromfile(filename, dtype=dt)
    df = pandas.DataFrame.from_records(data)
    print ('...Done')
    return np.array(df.T)[0]

data = read_data('qoi_M1024.dat')
#data = generate_celerite_data(100000)
n_max = min(len(data),100000)
tauint = []
n_chunk = 10
for j in range(n_chunk):
    autocorr = AutoCorrelation(data[j*n_max:(j+1)*n_max])
    tauint.append(autocorr.tau_int())

avg = np.average(tauint)
error = np.sqrt(1./(n_chunk-1.)*np.var(tauint))
print (avg,' +/- ',error)
autocorr = AutoCorrelation(data[4*n_max:5*n_max])
autocorr.plot('correlation.pdf')

from matplotlib import pyplot as plt
import numpy as np
import scipy.special
import math
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show

a = 0.1
I = 1.0

def f(phi_p,phi_m,phi):
    '''Function in exponential'''
    return I/a*(2.-math.cos(phi_p-phi)-math.cos(phi-phi_m))

def C(phi_p,phi_m):
    '''Constant C'''
    return 2.*math.fabs(math.cos(0.5*(phi_p-phi_m)))

def delta_phi(phi_p,phi_m):
    '''Shift delta phi'''
    return math.atan2(math.sin(phi_m)+math.sin(phi_p),math.cos(phi_m)+math.cos(phi_p))

def f_approx(phi_p,phi_m,phi):
    '''Approximate function f'''
    dphi = delta_phi(phi_p,phi_m)
    phi_2 = phi
    if (phi_2 > dphi+math.pi):
        phi_2 -= 2.*math.pi
    if (phi_2 < dphi-math.pi):
        phi_2 += 2.*math.pi
    return I/a*(2-C(phi_p,phi_m)+0.5*C(phi_p,phi_m)*(phi_2-dphi)**2)

def prob(g,norm,phi_p,phi_m):
    ''' Convert to probability by exponentiating function g'''
    dphi = delta_phi(phi_p,phi_m)
    nrm = norm(phi_p,phi_m)
    return lambda phi: nrm*math.exp(-g(phi_p,phi_m,phi)+g(phi_p,phi_m,dphi))

def Znorm(phi_p,phi_m):
    '''Normalisation factor Z^{-1}'''
    sigma = 2*I*C(phi_p,phi_m)/a
    return 1./(2.*math.pi*math.exp(-0.5*sigma)*scipy.special.iv(0,0.5*sigma))

def Znorm_approx(phi_p,phi_m):
    '''Approximate normalisation factor Z^{-1}'''
    sigma = 2*I*C(phi_p,phi_m)/a
    return math.sqrt(sigma/(4.*math.pi))

#phi_m = 0.0
#phi_p = 0.99*math.pi

phi_m = 0.2
phi_p = 1.8

print 'delta phi = ',delta_phi(phi_p,phi_m)
print 'C         = ',C(phi_p,phi_m)
print 'C^{-1/2}  = ',C(phi_p,phi_m)**(-0.5)

X = np.arange(-math.pi,math.pi,0.001)

plots = []
legends = []

plt.clf()
ax = plt.gca()
ax.set_xlim(-math.pi,math.pi)
ax.set_xlabel(r'$\phi$')
ax.set_xticks((-math.pi,-0.5*math.pi,0,0.5*math.pi,math.pi))
ax.set_xticklabels((r'$-\pi$',r'$-\frac{\pi}{2}$',r'$0$',r'$\frac{\pi}{2}$',r'$\pi$'))
plots.append(plt.plot(X,np.vectorize(lambda x: prob(f,Znorm,phi_p,phi_m)(x))(X),
                      linewidth=2,color='blue')[0])
legends.append(r'$\mathcal{Z}\cdot p(\phi|\phi_+,\phi_-)$')
plots.append(plt.plot(X,np.vectorize(lambda x: prob(f_approx,Znorm_approx,phi_p,phi_m)(x))(X),
                      linewidth=2,color='red')[0])
legends.append(r'$\mathcal{Z}\cdot p_{\operatorname{approx}}(\phi|\phi_+,\phi_-)$')
plots.append(plt.plot(X,np.vectorize(lambda x: prob(f_approx,Znorm_approx,phi_p,phi_m)(x)/prob(f,Znorm,phi_p,phi_m)(x))(X),
                      linewidth=2,color='green')[0])
legends.append(r'ratio')
plt.legend(plots,legends)
plt.savefig('rotor_probability.pdf',bbox_inches='tight')

plt.clf()
ax = plt.gca()
ax.set_xlim(-math.pi,math.pi)
ax.set_xlabel(r'$\phi$')
ax.set_xticks((-math.pi,-0.5*math.pi,0,0.5*math.pi,math.pi))
ax.set_xticklabels((r'$-\pi$',r'$-\frac{\pi}{2}$',r'$0$',r'$\frac{\pi}{2}$',r'$\pi$'))
plots = []
labels = []
plots.append(plt.plot(X,np.vectorize(lambda x: f(phi_p,phi_m,x))(X),
                      linewidth=2,color='blue')[0])
legends.append(r'$f(\phi;\phi_+,\phi_-)$')
plots.append(plt.plot(X,np.vectorize(lambda x: f_approx(phi_p,phi_m,x))(X),
                      linewidth=2,color='red')[0])
legends.append(r'$f_{\operatorname{approx}}(\phi;\phi_+,\phi_-)$')
plt.legend(plots,legends)
plt.savefig('rotor_log_probability.pdf',bbox_inches='tight')

import numpy as np
from scipy import linalg as la

class GFFAction(object):
    def __init__(self,Mlat,alpha,beta,h_offset=None):
        self.Mlat = Mlat
        self.Mlat_coarse = self.Mlat//2
        self.ndof = self.Mlat**2
        self.alpha = alpha
        self.beta = beta
        self._build_cholesky()
        self.h_offset = h_offset
        if not self.h_offset is None:
            self.set_h_offset(self.h_offset)

    def _cart2lin_idx(self,i,j):
        return (i % self.Mlat) + self.Mlat*(j % self.Mlat)

    def _lin2cart_idx(self,ell):
        i = ell % self.Mlat
        j = ell // self.Mlat
        return i,j

    def _build_cholesky(self):
        '''Construct Cholesky factorisation for sampling'''
        self.Q = np.zeros((self.ndof,self.ndof))
        for i in range(self.Mlat):
            for j in range(self.Mlat):
                ell =  self._cart2lin_idx(i,j)
                self.Q[ell,ell] = 4.*self.alpha+self.beta
                for offset in ((+1,0),(-1,0),(0,+1),(0,-1)):
                    ell_prime = self._cart2lin_idx(i+offset[0],j+offset[1])
                    self.Q[ell,ell_prime] = -self.alpha
        self.cholesky_LT = np.linalg.cholesky(self.Q).transpose()

    def set_h_offset(self,h_offset):
        self.h_offset = h_offset
        self.Q_inv_h_offset = np.linalg.solve(self.Q,h_offset.reshape(self.ndof)).reshape((self.Mlat,self.Mlat))

    def draw(self,phi_state):
        '''Draw the unknowns at all sites using the Cholesky factorisation of
        the precision matrix.

        :arg phi_state: State to populate
        '''
        psi = np.random.normal(size=self.ndof)
        phi_state[:,:] = la.solve_triangular(self.cholesky_LT,psi).reshape((self.Mlat,self.Mlat))
        if not self.h_offset is None:
            phi_state[:,:] -= self.Q_inv_h_offset[:,:]

    def gibbs_update(self,phi_state,ell):
        '''Gibbs update of a site with given linear index'''
        i,j = self._lin2cart_idx(ell)
        Delta = self.alpha * ( phi_state[(i+1)%self.Mlat,j%self.Mlat]
                             + phi_state[(i-1)%self.Mlat,j%self.Mlat]
                             + phi_state[i%self.Mlat,(j+1)%self.Mlat]
                             + phi_state[i%self.Mlat,(j-1)%self.Mlat] )
        if not self.h_offset is None:
            Delta -= self.h_offset[i,j]
        mu = Delta / (4.*self.alpha+self.beta)
        sigma = 1./np.sqrt(4.*self.alpha+self.beta)
        phi_state[i,j] = np.random.normal(loc=mu,scale=sigma)

    def restrict(self,phi_state,coarse_action=None):
        '''Construct coarsened version of action for conditional sampling'''
        h_coarse = np.zeros((self.Mlat_coarse,self.Mlat_coarse))
        for i in range(self.Mlat_coarse):
            for j in range(self.Mlat_coarse):
                h_ij = (2.*self.alpha+self.beta)*( phi_state[(2*i-1)%self.Mlat,(2*j-1)%self.Mlat]
                                              + phi_state[(2*i-1)%self.Mlat,(2*j+1)%self.Mlat]
                                              + phi_state[(2*i+1)%self.Mlat,(2*j-1)%self.Mlat]
                                              + phi_state[(2*i+1)%self.Mlat,(2*j+1)%self.Mlat] )
                h_ij -= self.alpha * ( phi_state[(2*i-3)%self.Mlat,(2*j-1)%self.Mlat]
                                         + phi_state[(2*i-3)%self.Mlat,(2*j+1)%self.Mlat]
                                         + phi_state[(2*i+3)%self.Mlat,(2*j-1)%self.Mlat]
                                         + phi_state[(2*i+3)%self.Mlat,(2*j+1)%self.Mlat]
                                         + phi_state[(2*i-1)%self.Mlat,(2*j-3)%self.Mlat]
                                         + phi_state[(2*i+1)%self.Mlat,(2*j-3)%self.Mlat]
                                         + phi_state[(2*i-1)%self.Mlat,(2*j+3)%self.Mlat]
                                         + phi_state[(2*i+1)%self.Mlat,(2*j+3)%self.Mlat] )
                if not self.h_offset is None:
                    h_ij +=  ( self.h_offset[(2*i-1)%self.Mlat,(2*j-1)%self.Mlat]
                             + self.h_offset[(2*i-1)%self.Mlat,(2*j+1)%self.Mlat]
                             + self.h_offset[(2*i+1)%self.Mlat,(2*j-1)%self.Mlat]
                             + self.h_offset[(2*i+1)%self.Mlat,(2*j+1)%self.Mlat] )
                h_coarse[i,j] = h_ij
        if coarse_action is None:
            return GFFAction(self.Mlat_coarse,2.*self.alpha,4.*self.beta,h_coarse)
        else:
            coarse_action.set_h_offset(h_coarse)

class TwogridMonteCarlo(object):
    def __init__(self,action,npresmooth,npostsmooth):
        self.action = action
        self.Mlat = self.action.Mlat
        self.ndof = self.action.ndof
        self.coarse_action = None
        self.npresmooth = npresmooth
        self.npostsmooth = npostsmooth
        self.dof_indices = np.arange(self.ndof)
        self.phi_coarse = np.zeros((self.Mlat//2,self.Mlat//2))

    def draw(self,phi_state):
        for k in range(self.npresmooth):
            self.gibbs_update(phi_state)
        if self.coarse_action is None:
            self.coarse_action = self.action.restrict(phi_state)
        else:
            self.action.restrict(phi_state,self.coarse_action)
        self.coarse_action.draw(self.phi_coarse)
        self.prolong_add(self.phi_coarse,phi_state)
        for k in range(self.npostsmooth):
            self.gibbs_update(phi_state)

    def prolong_add(self,phi_coarse,phi_fine):
        '''Prolongate from coarse to fine grid'''
        for i in range(self.coarse_action.Mlat):
            for j in range(self.coarse_action.Mlat):
                phi_coarse_ij = phi_coarse[i,j]
                phi_fine[2*i,  2*j  ] += phi_coarse_ij
                phi_fine[2*i+1,2*j  ] += phi_coarse_ij
                phi_fine[2*i  ,2*j+1] += phi_coarse_ij
                phi_fine[2*i+1,2*j+1] += phi_coarse_ij

    def gibbs_update(self,phi_state):
        '''Carry out a Gibbs sweep over the entire lattice'''
        np.random.shuffle(self.dof_indices)
        for ell in self.dof_indices:
            self.action.gibbs_update(phi_state,ell)


class QoISquaredField(object):
    def __init__(self,Mlat,mass):
        self.Mlat = Mlat
        self.mass = mass

    def evaluate(self,phi_state):
        return np.mean(phi_state[:,:]**2)
        S = 0.0
        for i in range(self.Mlat):
            for j in range(self.Mlat):
                S += phi_state[i,j]*phi_state[(i+self.Mlat//4)%self.Mlat,(j+self.Mlat//4)%self.Mlat]
        return S/(self.Mlat**2)


    def exact_value(self):
        S = 0
        mu2 = (self.mass/self.Mlat)**2
        for k1 in range(self.Mlat):
            for k2 in range(self.Mlat):
                S += 1./(4.*(np.sin(np.pi*k1/self.Mlat)**2+np.sin(np.pi*k2/self.Mlat)**2)+mu2)
        return S/self.Mlat**2

if __name__ == '__main__':
    Mlat = 16
    mass = 10.
    nsamples = 10000
    alat = 1./Mlat
    mu2 = (alat*mass)**2
    npresmooth = 2
    npostsmooth = 2
    action = GFFAction(Mlat,1.,mu2)
    tgmc = TwogridMonteCarlo(action,npresmooth,npostsmooth)
    qoi = QoISquaredField(Mlat,mass)
    qoi_coarse = QoISquaredField(Mlat//2,mass)


    # ==== Single-level method ====
    phi_state = np.zeros((Mlat,Mlat))
    S = 0
    S2 = 0
    S_coarse = 0.0
    S2_coarse = 0.0
    for k in range(nsamples):
        tgmc.draw(phi_state)
        q = qoi.evaluate(phi_state)
        S += q
        S2 += q**2
        q_coarse = qoi_coarse.evaluate(tgmc.phi_coarse)
        S_coarse += q_coarse
        S2_coarse += q_coarse**2
    average = S/nsamples
    variance = np.sqrt(1./(nsamples-1)*(S2/nsamples-(S/nsamples)**2))
    coarse_average = S_coarse/nsamples
    coarse_variance = np.sqrt(1./(nsamples-1)*(S2_coarse/nsamples-(S_coarse/nsamples)**2))
    print ('numerical [fine]   = ', ('%8.4f' % average),' +/- ',('%8.4f' % variance))
    print ('numerical [coarse] = ', ('%8.4f' % coarse_average),' +/- ',('%8.4f' % coarse_variance))
    print ('analytical         = ', ('%8.4f' % qoi.exact_value()))

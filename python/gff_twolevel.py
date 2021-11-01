import numpy as np
from scipy import linalg as la

class GFFAction(object):
    def __init__(self,Mlat,mass):
        self.exact_schur = False
        self.Mlat = Mlat
        self.alat = 1./self.Mlat
        self.mass = mass
        self.ndof = self.Mlat**2
        self._build_cholesky()
        self._build_coarse_cholesky()

    def _cart2lin_idx(self,i,j):
        return (i % self.Mlat) + self.Mlat*(j % self.Mlat)

    def _build_cholesky(self):
        '''Construct (upper triangular) Cholesky factor L^T of the precision
        matrix Q on the fine level.
        '''
        Q = np.zeros((self.ndof,self.ndof))
        # coarse level non-dimensionaled squared mass
        mu2 = (self.alat*self.mass)**2
        for i in range(self.Mlat):
            for j in range(self.Mlat):
                ell =  self._cart2lin_idx(i,j)
                Q[ell,ell] = 4.+mu2
                for offset in ((+1,0),(-1,0),(0,+1),(0,-1)):
                    ell_prime = self._cart2lin_idx(i+offset[0],j+offset[1])
                    Q[ell,ell_prime] = -1.0
        self.cholesky_LT = np.linalg.cholesky(Q).transpose()

    def _build_coarse_cholesky(self):
        '''Construct (upper triangular) Cholesky factor L^T of the precision
        matrix Q used for coarse level sampling. Note that the matrix Q is
        defined on *all* lattice sites, but it is the identity matrix on the
        non-coarse sites.
        '''
        Q = np.zeros((self.ndof,self.ndof))
        # coarse level non-dimensionaled squared mass
        mu2_coarse = 2.*(self.alat*self.mass)**2
        for i in range(self.Mlat):
            for j in range(self.Mlat):
                ell =  self._cart2lin_idx(i,j)
                if (i+j)%2==0:
                    Q[ell,ell] = 4.+mu2_coarse
                    for offset in ((+1,+1),(+1,-1),(-1,+1),(-1,-1)):
                        ell_prime = self._cart2lin_idx(i+offset[0],j+offset[1])
                        Q[ell,ell_prime] = -1.0
                else:
                    Q[ell,ell] = 1.0
        self.cholesky_coarse_LT = np.linalg.cholesky(Q).transpose()

    def evaluate_coarse(self,phi_state):
        '''Evaluate the coarse-level action at all the coarse sites only

        :arg phi_state: state to evaluate
        '''
        S = 0.0
        mu2_coarse = 2.*(self.alat*self.mass)**2
        for i in range(self.Mlat):
            for j in range(self.Mlat):
                if (i+j)%2==0:
                    phi_local = phi_state[i,j]
                    S_local = (4.+mu2_coarse)*phi_local
                    for offset in ((+1,+1),(+1,-1),(-1,+1),(-1,-1)):
                        S_local -= phi_state[(i+offset[0])%self.Mlat,(j+offset[1])%self.Mlat]
                    S += S_local*phi_local
        return 0.5*S

    def evaluate(self,phi_state):
        '''Evaluate the fine-level action at all sites

        :arg phi_state: state to evaluate
        '''
        S = 0.0
        mu2 = (self.alat*self.mass)**2
        for i in range(self.Mlat):
            for j in range(self.Mlat):
                phi_local = phi_state[i,j]
                S_local = (4.+mu2)*phi_local
                for offset in ((+1,0),(-1,0),(0,+1),(0,-1)):
                    S_local -= phi_state[(i+offset[0])%self.Mlat,(j+offset[1])%self.Mlat]
                S += S_local*phi_local
        return 0.5*S

    def evaluate_fillin(self,phi_state):
        '''Evaluate the total value of the fill-in action at the edges

        :arg phi_state: state to evaluate
        '''
        sigma2_inv = 4.+(self.alat*self.mass)**2
        sigma2 = 1./sigma2_inv
        S = 0.0
        for i in range(self.Mlat):
            for j in range(self.Mlat):
                if (i+j)%2==1:
                    Delta = 0.0
                    for offset in ((+1,0),(-1,0),(0,+1),(0,-1)):
                        Delta += phi_state[(i+offset[0])%self.Mlat,(j+offset[1])%self.Mlat]
                    S += (phi_state[i,j]-sigma2*Delta)**2
        return 0.5*sigma2_inv*S

    def draw(self,phi_state):
        '''Draw the unknowns at all sites using the Cholesky factorisation of
        the precision matrix.

        :arg phi_state: State to populate
        '''
        psi = np.random.normal(size=self.ndof)
        phi_state[:,:] = la.solve_triangular(self.cholesky_LT,psi).reshape((self.Mlat,self.Mlat))

    def draw_coarse_dofs(self,phi_state):
        '''Draw the unknowns at the coarse lattice sites from the distribution
        defined by the coarse action.

        :arg phi_state: State to populate
        '''
        psi = np.random.normal(size=self.ndof)
        phi_state[:,:] = la.solve_triangular(self.cholesky_coarse_LT,psi).reshape(self.Mlat,self.Mlat)

    def draw_fillin_dofs(self,phi_state):
        '''Draw the unknowns at the edge lattice sites from the fill-in
        distribution

        :arg phi_state: State to populate
        '''
        sigma2 = 1./(4.+(self.alat*self.mass)**2)
        sigma = np.sqrt(sigma2)
        for i in range(self.Mlat):
            for j in range(self.Mlat):
                if (i+j)%2==1:
                    Delta = 0.0
                    for offset in ((+1,0),(-1,0),(0,+1),(0,-1)):
                        Delta += phi_state[(i+offset[0])%self.Mlat,(j+offset[1])%self.Mlat]
                    phi_state[i,j] = np.random.normal(sigma2*Delta,sigma)

class QoISquaredField(object):
    def __init__(self,Mlat,mass,coarse_only=False):
        self.Mlat = Mlat
        self.mass = mass
        self.coarse_only = coarse_only

    def evaluate(self,phi_state):
        if (self.coarse_only):
            q = 0
            for i in range(self.Mlat/2):
                for j in range(self.Mlat/2):
                    q += phi_state[2*i,2*j]**2
            return 4.*q/(self.Mlat**2)
        else:
            return np.mean(phi_state**2)

    def exact_value(self):
        S = 0
        mu2 = (self.mass/self.Mlat)**2
        for k1 in range(self.Mlat):
            for k2 in range(self.Mlat):
                S += 1./(4.*(np.sin(np.pi*k1/self.Mlat)**2+np.sin(np.pi*k2/self.Mlat)**2)+mu2)
        return S/self.Mlat**2

class TwoLevelSampler(object):
    def __init__(self,action,qoi):
        self.action = action
        self.qoi = qoi
        self.Mlat = self.action.Mlat
        self.phi_current = np.zeros((self.Mlat,self.Mlat))
        self.phi_proposal = np.zeros((self.Mlat,self.Mlat))
        self.action.draw(self.phi_current)
        self.S_fine = self.action.evaluate(self.phi_current)
        self.S_coarse = self.action.evaluate_coarse(self.phi_current)
        self.S_fillin = self.action.evaluate_fillin(self.phi_current)
        self.nsamples = 0
        self.naccepted = 0

    def step(self):
        self.action.draw_coarse_dofs(self.phi_proposal)
        self.action.draw_fillin_dofs(self.phi_proposal)
        S_fine_proposal = self.action.evaluate(self.phi_proposal)
        S_coarse_proposal = self.action.evaluate_coarse(self.phi_proposal)
        S_fillin_proposal = self.action.evaluate_fillin(self.phi_proposal)
        DeltaS = 0
        DeltaS += S_fine_proposal - self.S_fine
        DeltaS += self.S_coarse - S_coarse_proposal
        DeltaS += self.S_fillin - S_fillin_proposal
        accept = False
        if (DeltaS < 0):
            accept = True
        else:
            accept = (np.random.uniform(0.0,1.0) <= np.exp(-DeltaS))
        if accept:
            self.naccepted += 1
            self.S_fine = S_fine_proposal
            self.S_coarse = S_coarse_proposal
            self.S_fillin = S_fillin_proposal
        self.nsamples += 1

    def rho_accept(self):
        return self.naccepted/self.nsamples

Mlat = 8
mass = 10.0
nsamples = 10000
action = GFFAction(Mlat,mass)
qoi = QoISquaredField(Mlat,mass)

# ==== Single-level method ====
phi_state = np.zeros((Mlat,Mlat))
S = 0
for k in range(nsamples):
    action.draw(phi_state)
    S += qoi.evaluate(phi_state)
print ('numerical  = ', S/nsamples)
print ('analytical = ', qoi.exact_value())

# ==== Two-level method ===
sampler = TwoLevelSampler(action,qoi)
for k in range(nsamples):
    sampler.step()
print ('total samples    = ',sampler.nsamples)
print ('accepted samples = ',sampler.naccepted)
print ('acceptance rate  = ',sampler.rho_accept())
print ('rejection rate   = ',1.-sampler.rho_accept())

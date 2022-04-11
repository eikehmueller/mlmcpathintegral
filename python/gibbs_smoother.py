class CoarseGibbsSmoother(object):
    def __init__(self,action,nsmooth=1):
        self.Mlat = action.Mlat
        self.alat = action.alat
        self.ndof = self.Mlat**2//2
        mu2 = (action.mass*action.alat)**2
        self.stencil_initial = {(0,0): 4.+2.*mu2,
                                (+1,+1):-1.,
                                (+1,-1):-1.,
                                (-1,+1):-1.,
                                (-1,-1):-1.}
        tmp = 1./(4.+mu2)
        self.stencil = {(0,0):4.+mu2-4.*tmp,
                        (+1,+1):-2*tmp,
                        (+1,-1):-2*tmp,
                        (-1,+1):-2*tmp,
                        (-1,-1):-2*tmp,
                        (+2,0):-tmp,
                        (-2,0):-tmp,
                        (0,+2):-tmp,
                        (0,-2):-tmp}
        self.nsmooth = nsmooth
        self._build_index_maps()
        self._build_matrices()

    def _build_index_maps(self):
        ell = 0
        self.cart2lin_idx = {}
        self.lin2cart_idx = {}
        for i in range(self.Mlat):
            for j in range(self.Mlat):
                if (i+j)%2==0:
                    self.lin2cart_idx[ell] = (i,j)
                    self.cart2lin_idx[(i,j)] = ell
                    ell+=1

    def _build_matrices(self):
        Q_prec = np.zeros((self.ndof,self.ndof))
        Q_prec_initial = np.zeros((self.ndof,self.ndof))
        self.Sigma_initial = np.zeros((self.ndof,self.ndof))
        self.Sigma = np.zeros((self.ndof,self.ndof))
        # build precision matrices
        for ell in range(self.ndof):
            (i,j) = self.lin2cart_idx[ell]
            for offset in self.stencil_initial.keys():
                ell_offset = self.cart2lin_idx[((i+offset[0])%self.Mlat,(j+offset[1])%self.Mlat)]
                Q_prec_initial[ell,ell_offset] = self.stencil_initial[offset]
            for offset in self.stencil.keys():
                ell_offset = self.cart2lin_idx[((i+offset[0])%self.Mlat,(j+offset[1])%self.Mlat)]
                Q_prec[ell,ell_offset] = self.stencil[offset]
        # resulting covariance matrices
        self.Sigma = np.linalg.inv(Q_prec)
        self.Sigma_initial = np.linalg.inv(Q_prec_initial)
        M = np.tril(Q_prec,k=0)
        G = np.identity(self.ndof) - np.linalg.inv(M)@Q_prec
        G_k = np.linalg.matrix_power(G,self.nsmooth)
        self.Sigma_iter = self.Sigma + G_k@(self.Sigma_initial-self.Sigma)@G_k.transpose()
        self.Q_prec_iter = np.linalg.inv(self.Sigma_iter)

    def apply(self,phi_state):
        a_diag = self.stencil[(0,0)]
        inv_sqrt_a_diag = 1./np.sqrt(a_diag)
        stencil_offdiag = dict(self.stencil)
        stencil_offdiag.pop((0,0))
        for k in range(self.nsmooth):
            for i in range(self.Mlat):
                for j in range(self.Mlat):
                    if (i+j)%2==0:
                        mu = 0.0
                        for offset in stencil_offdiag:
                            mu += stencil_offdiag[offset]*phi_state[(i+offset[0])%self.Mlat,
                                                                    (j+offset[1])%self.Mlat]
                        mu /= a_diag
                        phi_state[i,j] = np.random.normal(loc=mu, scale=inv_sqrt_a_diag)

    def evaluate(self,phi_state):
        phi_lin = np.zeros(self.ndof)
        for i in range(self.Mlat):
            for j in range(self.Mlat):
                if (i+j)%2==0:
                    ell = self.cart2lin_idx[(i,j)]
                    phi_lin[ell] = phi_state[i,j]
        return 0.5*phi_lin.transpose()@self.Q_prec_iter@phi_lin

    def _cov_function(self,cov):
        X = np.zeros(self.Mlat)
        Y = np.zeros(self.Mlat)
        for j in range(self.Mlat):
            ell = self.cart2lin_idx[(j,j)]
            X[j] = np.sqrt(2.)*self.alat*j
            Y[j] = cov[0,ell]
        return X,Y

    def covariance_function_initial(self):
        return self._cov_function(self.Sigma_initial)

    def covariance_function_iter(self):
        return self._cov_function(self.Sigma_iter)

    def covariance_function(self):
        return self._cov_function(self.Sigma)
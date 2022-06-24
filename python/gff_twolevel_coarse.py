import numpy as np
from matplotlib import pyplot as plt
from scipy import linalg as la
from gibbs_smoother import *


class GFFAction(object):
    """Action for Gaussian free field"""

    def __init__(self, Mlat, mass):
        self.exact_schur = False
        self.Mlat = Mlat
        self.alat = 1.0 / self.Mlat
        self.mass = mass
        self.ndof = self.Mlat**2
        self._build_cholesky()
        self._build_projector()

    def _cart2lin_idx(self, i, j):
        return (i % self.Mlat) + self.Mlat * (j % self.Mlat)

    def _build_cholesky(self):
        """Construct (upper triangular) Cholesky factor L^T of the precision
        matrix Q.
        """
        # 5-point stencil
        Q_5pt = np.zeros((self.ndof, self.ndof))
        mu2 = (self.alat * self.mass) ** 2
        for i in range(self.Mlat):
            for j in range(self.Mlat):
                ell = self._cart2lin_idx(i, j)
                Q_5pt[ell, ell] = 4.0 + 2.0 * mu2
                for offset in ((+1, 0), (-1, 0), (0, +1), (0, -1)):
                    ell_prime = self._cart2lin_idx(i + offset[0], j + offset[1])
                    Q_5pt[ell, ell_prime] = -1.0
        self.cholesky_5pt_LT = np.linalg.cholesky(Q_5pt).transpose()
        # 9-point stencil
        C_diag = (4.0 + mu2) - 4.0 / (4.0 + mu2)  # diagonal term (self-coupling)
        C_nn = -2.0 / (4.0 + mu2)  # nearest neighbour coupling
        C_nnn = -1.0 / (4.0 + mu2)  # next-to-nearest neighbour coupling
        Q_9pt = np.zeros((self.ndof, self.ndof))
        mu2 = (self.alat * self.mass) ** 2
        for i in range(self.Mlat):
            for j in range(self.Mlat):
                ell = self._cart2lin_idx(i, j)
                Q_9pt[ell, ell] = C_diag
                for offset in ((+1, 0), (-1, 0), (0, +1), (0, -1)):
                    ell_prime = self._cart2lin_idx(i + offset[0], j + offset[1])
                    Q_9pt[ell, ell_prime] = C_nn
                for offset in ((+1, +1), (+1, -1), (-1, +1), (-1, -1)):
                    ell_prime = self._cart2lin_idx(i + offset[0], j + offset[1])
                    Q_9pt[ell, ell_prime] = C_nnn
        self.cholesky_9pt_LT = np.linalg.cholesky(Q_9pt).transpose()

    def _build_projector(self):
        """Construct projection matrix"""
        U = np.zeros((self.Mlat, self.Mlat))
        for n in range(self.Mlat):
            for j in range(self.Mlat):
                if n == 0:
                    U[n, j] = 1 / np.sqrt(self.Mlat)
                else:
                    U[n, j] = np.sqrt(2 / self.Mlat) * np.cos(
                        np.pi / self.Mlat * (j + 0.5) * n
                    )
        D = np.zeros((self.Mlat, self.Mlat))
        for j in range(self.Mlat // 4):
            D[j, j] = 1.0
        self.P = U.T @ D @ U

    def project(self, phi_state, phi_projected):
        """project a field

        :arg phi_state: input field
        :arg phi_projected: projected output field"""
        phi_projected[:, :] = self.P @ phi_state[:, :] @ self.P.T

    def evaluate_5pt(self, phi_state):
        """Evaluate the action with a 5-point stencil

        :arg phi_state: state to evaluate
        """
        S = 0.0
        mu2 = (self.alat * self.mass) ** 2
        for i in range(self.Mlat):
            for j in range(self.Mlat):
                phi_local = phi_state[i, j]
                S_local = (4.0 + 2.0 * mu2) * phi_local
                for offset in ((+1, 0), (-1, 0), (0, +1), (0, -1)):
                    S_local -= phi_state[
                        (i + offset[0]) % self.Mlat, (j + offset[1]) % self.Mlat
                    ]
                S += S_local * phi_local
        return 0.5 * S

    def evaluate_9pt(self, phi_state):
        """Evaluate the action with a 9-point stencil

        :arg phi_state: state to evaluate
        """
        S = 0.0
        mu2 = (self.alat * self.mass) ** 2
        C_diag = (4.0 + mu2) - 4.0 / (4.0 + mu2)  # diagonal term (self-coupling)
        C_nn = -2.0 / (4.0 + mu2)  # nearest neighbour coupling
        C_nnn = -1.0 / (4.0 + mu2)  # next-to-nearest neighbour coupling
        for i in range(self.Mlat):
            for j in range(self.Mlat):
                phi_local = phi_state[i, j]
                S_local = C_diag * phi_local
                for offset in ((+1, 0), (-1, 0), (0, +1), (0, -1)):
                    S_local += (
                        C_nn
                        * phi_state[
                            (i + offset[0]) % self.Mlat, (j + offset[1]) % self.Mlat
                        ]
                    )
                for offset in ((+1, +1), (+1, -1), (-1, +1), (-1, -1)):
                    S_local += (
                        C_nnn
                        * phi_state[
                            (i + offset[0]) % self.Mlat, (j + offset[1]) % self.Mlat
                        ]
                    )
                S += S_local * phi_local
        return 0.5 * S

    def draw_5pt(self, phi_state):
        """Draw the unknowns at all sites using the Cholesky factorisation of
        the precision matrix for the 5 point stencil.

        :arg phi_state: State to populate
        """
        psi = np.random.normal(size=self.ndof)
        phi_state[:, :] = la.solve_triangular(self.cholesky_5pt_LT, psi).reshape(
            (self.Mlat, self.Mlat)
        )

    def draw_9pt(self, phi_state):
        """Draw the unknowns at all sites using the Cholesky factorisation of
        the precision matrix for the 9 point stencil.

        :arg phi_state: State to populate
        """
        psi = np.random.normal(size=self.ndof)
        phi_state[:, :] = la.solve_triangular(self.cholesky_9pt_LT, psi).reshape(
            (self.Mlat, self.Mlat)
        )


class TwoLevelSampler(object):
    def __init__(self, action):
        self.action = action
        self.nsamples = 0
        self.naccepted = 0
        self.phi_state = np.zeros((self.action.Mlat, self.action.Mlat))
        self.phi_proposal = np.zeros((self.action.Mlat, self.action.Mlat))
        self.action.draw_9pt(self.phi_state)

    def step(self):
        self.action.draw_5pt(self.phi_proposal)
        self.action.draw_9pt(self.phi_state)
        project = True
        if project:
            phi_state_projected = np.zeros((self.action.Mlat, self.action.Mlat))
            phi_proposal_projected = np.zeros((self.action.Mlat, self.action.Mlat))
            self.action.project(self.phi_state, phi_state_projected)
            self.action.project(self.phi_proposal, phi_proposal_projected)
        else:
            phi_state_projected = self.phi_state
            phi_proposal_projected = self.phi_proposal
        deltaS = (
            self.action.evaluate_9pt(phi_proposal_projected)
            - self.action.evaluate_9pt(phi_state_projected)
            + self.action.evaluate_5pt(phi_state_projected)
            - self.action.evaluate_5pt(phi_proposal_projected)
        )
        if deltaS < 0:
            accept = True
        else:
            accept = np.random.uniform(low=0.0, high=1.0) < np.exp(-deltaS)
        if accept:
            self.phi_state[:, :] = self.phi_proposal[:, :]
            self.naccepted += 1
        self.nsamples += 1

    @property
    def rho_accept(self):
        return self.naccepted / self.nsamples


Mlat = 64
mass = 1.0
nsamples = 100
action = GFFAction(Mlat, mass)

# ==== Two-level method ===
sampler = TwoLevelSampler(action)
for k in range(nsamples):
    sampler.step()
print("total samples    = ", sampler.nsamples)
print("accepted samples = ", sampler.naccepted)
print("acceptance rate  = {:4.2f}%".format(100 * sampler.rho_accept))
print("rejection rate   = {:4.2f}%".format(100 * (1.0 - sampler.rho_accept)))

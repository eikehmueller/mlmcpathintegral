import numpy as np
from matplotlib import pyplot as plt
from scipy import linalg as la
from gibbs_smoother import *


class Lattice2d(object):
    """A structured two-dimensional lattice"""

    def __init__(self, Mlatx, Mlaty, fine_lattice=None):
        """Construct new instance

        :arg Mlatx: Number of lattice points in x-direction
        :arg Mlaty: Number of lattice points in y-direction
        :arg fine_lattice: Next-finer lattice in the hierarchy
        """
        self.Mlatx = Mlatx
        self.Mlaty = Mlaty
        self.alatx = 1 / self.Mlatx
        self.alaty = 1 / self.Mlaty
        self.fine_lattice = fine_lattice

    def _cart2lin_idx(self, i, j):
        """Convert index pair to linear index

        :arg i: index in x-direction
        :arg j: index in y-direction
        """
        return (i % self.Mlatx) + self.Mlatx * (j % self.Mlaty)

    def copy2fine(self, phi_state, phi_fine_state):
        """Copy a field on the current lattice to the next finer lattice"""
        # Case 1: lattice has been refined in y-direction
        if (self.Mlatx == self.fine_lattice.Mlatx) and (
            2 * self.Mlaty == self.fine_lattice.Mlaty
        ):
            # Case 1: lattice has been refined in y-direction
            phi_fine_state[:, ::2] = phi_state[:, :]
        elif (2 * self.Mlatx == self.fine_lattice.Mlatx) and (
            self.Mlaty == self.fine_lattice.Mlaty
        ):
            # Case 2: lattice has been refined in x-direction
            phi_fine_state[::2, :] = phi_state[:, :]
        else:
            raise NotImplementedError("Can not copy field to fine lattice.")


class GFFAction(object):
    """Gaussian free field action"""

    def __init__(self, lattice):
        """Construct new instance

        :arg lattce: Underlying two dimensional lattice
        """
        self.lattice = lattice
        self.ndof = self.lattice.Mlatx * self.lattice.Mlaty

    def _build_cholesky(self):
        """Construct (upper triangular) Cholesky factor L^T of the precision
        matrix Q.
        """
        Q = np.zeros((self.ndof, self.ndof))
        # coarse level non-dimensionaled squared mass
        coeff_diag = self.alatx * self.alaty * self.mass**2 + 2.0 * (
            self.alatx / self.alaty + self.alaty / self.alatx
        )
        rho_x = self.alaty / self.alatx
        rho_y = self.alatx / self.alaty
        for i in range(self.Mlatx):
            for j in range(self.Mlaty):
                ell = self._cart2lin_idx(i, j)
                Q[ell, ell] = coeff_diag
                for offset in (+1, -1):
                    ell_prime = self._cart2lin_idx(i + offset, j)
                    Q[ell, ell_prime] = -rho_x
                    ell_prime = self._cart2lin_idx(i, j + offset)
                    Q[ell, ell_prime] = -rho_y
        self.cholesky_LT = np.linalg.cholesky(Q).transpose()

        # Construct row-wise precision matrices and Cholesky factorisations
        coeff_diag = self.alatx * self.alaty * self.mass**2 + 2.0 * (
            self.alatx / self.alaty + self.alaty / self.alatx
        )
        rho_x = self.alaty / self.alatx
        rho_y = self.alatx / self.alaty
        # x-direction
        Qrow_x = np.zeros((self.Mlatx, self.Mlatx))
        for i in range(self.Mlatx):
            Qrow_x[i, i] = coeff_diag
            for offset in (+1, -1):
                Qrow_x[i, (i + 1) % self.Mlatx] = -rho_x
                Qrow_x[i, (i - 1) % self.Mlatx] = -rho_x
        # y-direction
        Qrow_y = np.zeros((self.Mlaty, self.Mlaty))
        for j in range(self.Mlaty):
            Qrow_y[j, j] = coeff_diag
            for offset in (+1, -1):
                Qrow_y[j, (j + 1) % self.Mlaty] = -rho_y
                Qrow_y[j, (j - 1) % self.Mlaty] = -rho_y
        self.cholesky_Lrowx = np.linalg.cholesky(Qrow_x)
        self.cholesky_Lrowy = np.linalg.cholesky(Qrow_y)

    def evaluate(self, phi_state):
        """Evaluate the action at all sites

        :arg phi_state: state to evaluate
        """
        S = 0.0
        coeff_diag = self.alatx * self.alaty * self.mass**2 + 2.0 * (
            self.alatx / self.alaty + self.alaty / self.alatx
        )
        rho_x = self.alaty / self.alatx
        rho_y = self.alatx / self.alaty
        for i in range(self.Mlatx):
            for j in range(self.Mlaty):
                phi_local = phi_state[i, j]
                S_local = coeff_diag * phi_local
                for offset in (+1, -1):
                    S_local -= rho_x * phi_state[(i + offset) % self.Mlatx, j]
                    S_local -= rho_y * phi_state[i, (j + offset) % self.Mlaty]
                S += S_local * phi_local
        return 0.5 * S

    def draw_odd_lines(self, phi_state, direction="x"):
        """Draw all unknowns along lines with odd indices

        If direction is 'x', draw the unknowns located at indices (i,2*j+1), given
        the unknowns at indices (i,2*j) for i=0,1,...,Mlatx, j=0,1,...,Mlaty/2.
        If direction is 'y;, draw the unknowns located at indices (2*i+1,j), given
        the unknowns at indices (2*i,j) for i=0,1,...,Mlatx/2, j=0,1,...,Mlaty.

        :arg phi_state: Field to populate
        :arg direction: Direction of lines to fill, has to be 'x' or 'y'
        """
        assert (direction == "x") or (direction == "y")
        if direction == "x":
            rho_y = self.alatx / self.alaty
            psi_r = np.zeros(self.Mlatx)
            for j in range(1, self.Mlaty, 2):
                psi_r[:] = (
                    0.5
                    * rho_y
                    * (
                        phi_state[:, (j + 1) % self.Mlaty]
                        + phi_state[:, (j - 1) % self.Mlaty]
                    )
                )
                xi = np.random.normal(size=self.Mlatx) + la.solve_triangular(
                    self.cholesky_LrowxT, psi_r, lower=True
                )
                phi_state[:, j] = la.solve_triangular(
                    self.cholesky_LrowxT.transpose(), xi
                )
        else:
            rho_x = self.alaty / self.alatx
            psi_r = np.zeros(self.Mlaty)
            for i in range(1, self.Mlatx, 2):
                psi_r[:] = (
                    0.5
                    * rho_x
                    * (
                        phi_state[(i + 1) % self.Mlatx, :]
                        + phi_state[(i - 1) % self.Mlatx, :]
                    )
                )
                xi = np.random.normal(size=self.Mlaty) + la.solve_triangular(
                    self.cholesky_LrowyT, psi_r, lower=True
                )
                phi_state[i, :] = la.solve_triangular(
                    self.cholesky_LrowyT.transpose(), xi
                )

    def draw(self, phi_state):
        """Draw the unknowns at all sites using the Cholesky factorisation of
        the precision matrix.

        :arg phi_state: State to populate
        """
        psi = np.random.normal(size=self.ndof)
        phi_state[:, :] = la.solve_triangular(self.cholesky_LT, psi).reshape(
            (self.Mlatx, self.Mlaty)
        )

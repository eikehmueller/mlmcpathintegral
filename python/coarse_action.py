from itertools import product
import numpy as np
from sympy import Indexed, Add, Mul, Pow, expand


class Action:
    """Class representing an action"""

    def __init__(self, expr, phi, sidxs):
        self.phi = phi
        self.sidxs = sidxs
        self.dim = len(sidxs)
        self.expr = self.normal_shift(expr)

    def _min_phi_idx(self, expr, phi):
        """Parse the given expression recursively to find the offsets
        alpha_0,...,alpha_{k-1} with the smallest absolute values of
        all terms of the form phi[i_0+alpha_0,...,i_{k-1}+alpha_{k-1}].

        Returns [alpha_0,...,alpha_{k-1}] as a list

        :arg expr: expression to parse
        """
        cls = expr.func
        if cls is Add or cls is Mul or cls is Pow:
            # for Add, Mul or Pow return minimum over all subexpressions
            return np.asarray([self._min_phi_idx(arg, phi) for arg in expr.args]).min(
                axis=0
            )
        elif cls is Indexed:
            if expr.base == phi:
                # return minimum offset [alpha_0,...,alpha_{k-1}] for expressions of
                # the form phi[i_0+alpha_0,...,i_{k-1}+alpha_{k-1}]
                idxs = expr.indices
                return [idx.subs([(ell, 0) for ell in self.sidxs]) for idx in idxs]
        # ignore all other terms
        return self.dim * [np.infty]

    def _shift_index(self, expr, shift, phi):
        """Parse the given expression recursively and make the following replacements:

            phi[i_0+alpha_0,...,i_{k-1}+alpha_{k-1}]
            -> phi[i_0+alpha_0-s_0,...,i_{k-1}+alpha_{k-1}-s_{k-1}]

        for the coarse fields phi and

            S[i_0+alpha_0,...,i_{k-1}+alpha_{k-1}]
            -> S[i_0+alpha_0-2*s_0,...,i_{k-1}+alpha_{k-1}-2*s_{k-1}]

        for all other indexed symbols S.

        :arg expr: expression to parse
        :arg shift: shift vector [s_0,s_1,...,s_{k-1}]
        :arg phi: symbol phi
        """
        cls = expr.func
        if cls is Add or cls is Mul or cls is Pow:
            # for Add, Mul and Pow nodes, recursively parse subexpressions
            return cls(*[self._shift_index(arg, shift, phi) for arg in expr.args])
        elif cls is Indexed:
            return expr.base[[idx - sidx for idx, sidx in zip(expr.indices, shift)]]
        # return unchanged expression by default
        return expr

    def normal_shift(self, expr):
        """Shift such that lowest idx on phi is (0,0,...,0)"""
        nexpr = 0
        for arg in expand(expr).args:
            shift = self._min_phi_idx(arg, self.phi)
            nexpr += self._shift_index(arg, shift, self.phi)
        return nexpr


class Coarsener:
    """Class for coarsening a given expression"""

    def __init__(self, phi_fine, phi_coarse, psi, sidxs):
        self.phi_fine = phi_fine
        self.phi_coarse = phi_coarse
        self.psi = psi
        self.sidxs = sidxs
        self.dim = len(sidxs)

    def _min_phi_idx(self, expr):
        """Parse the given expression recursively to find the smallest offsets
        alpha_0,...,alpha_{k-1} of all terms of the form
        phi_coarse[i_0+alpha_0,...,i_{k-1}+alpha_{k-1}].

        Returns [alpha_0,...,alpha_{k-1}] as a list

        :arg expr: expression to parse
        """
        cls = expr.func
        if cls is Add or cls is Mul or cls is Pow:
            # for Add, Mul or Pow return minimum over all subexpressions
            return np.asarray([self._min_phi_idx(arg) for arg in expr.args]).min(axis=0)
        elif cls is Indexed:
            if expr.base == self.phi_coarse:
                # return minimum offset [alpha_0,...,alpha_{k-1}] for expressions of
                # the form phi[i_0+alpha_0,...,i_{k-1}+alpha_{k-1}]
                idxs = expr.indices
                return [idx.subs([(ell, 0) for ell in self.sidxs]) for idx in idxs]
        # ignore all other terms
        return self.dim * [np.infty]

    def _shift_index(self, expr, shift):
        """Parse the given expression recursively and make the following replacements:

            phi_coarse[i_0+alpha_0,...,i_{k-1}+alpha_{k-1}]
            -> phi_coarse[i_0+alpha_0-s_0,...,i_{k-1}+alpha_{k-1}-s_{k-1}]

        for the coarse fields phi_coarse and

            S[i_0+alpha_0,...,i_{k-1}+alpha_{k-1}]
            -> S[i_0+alpha_0-2*s_0,...,i_{k-1}+alpha_{k-1}-2*s_{k-1}]

        for all other indexed symbols S.

        :arg expr: expression to parse
        :arg phi: symbol phi
        :arg shift: shift vector [s_0,s_1,...,s_{k-1}]
        """
        cls = expr.func
        if cls is Add or cls is Mul or cls is Pow:
            # for Add, Mul and Pow nodes, recursively parse subexpressions
            return cls(*[self._shift_index(arg, shift) for arg in expr.args])
        elif cls is Indexed:
            offset = 1 if expr.base == self.phi_coarse else 2
            return expr.base[
                [idx - offset * sidx for idx, sidx in zip(expr.indices, shift)]
            ]
        # return unchanged expression by default
        return expr

    def _coarsen_idx(self, expr, offsets):
        """
        In all indexed symbols S[f_0(i_0),...,f_{d-1}(i_{d-1})] in the
        expression recursively replace

            f_j(i_j) -> f_j(2*i_j+delta_j)

        :arg expr: expression to parse
        :arg offsets: offsets delta_j
        """
        cls = expr.func
        if cls is Add or cls is Mul or cls is Pow:
            # for Add, Mul or Pow return minimum over all subexpressions
            return cls(*[self._coarsen_idx(arg, offsets) for arg in expr.args])
        elif cls is Indexed:
            idxs = expr.indices
            return expr.base[
                [
                    idx.subs(
                        [
                            (ell, 2 * ell + offset)
                            for ell, offset in zip(self.sidxs, offsets)
                        ]
                    )
                    for idx in idxs
                ]
            ]
        # leave all other terms untransformed
        return expr

    def _restrict(self, expr):
        """In the given expression, recursively replace:

            phi_fine[j_0,...,j_{d-1}]
            -> phi_coarse[j_0//2,...,j_{d-1}//2] + psi[j_0,...,j_{d-1}]

        :arg expr: expression to parse
        """
        cls = expr.func
        if cls is Add or cls is Mul or cls is Pow:
            # for Add, Mul or Pow return minimum over all subexpressions
            return cls(*[self._restrict(arg) for arg in expr.args])
        elif cls is Indexed:
            if expr.base == self.phi_fine:
                idxs = expr.indices
                return self.phi_coarse[[idx // 2 for idx in idxs]] + self.psi[idxs]
        # leave all other terms untransformed
        return expr

    def apply(self, expr):
        """
        Coarsen the given expression by applying two transformations:

        1. Split the (implicit) summation in each dimension into even and odd parts
        2. Replace the fine level fields by the coarse level fields plus a correction, i.e.
            phi_fine -> phi_coarse + psi

        :arg expr: expression to parse
        """

        # split summation into even/odd part and restrict fields
        coarse_expr = 0
        for offsets in product([0, 1], repeat=self.dim):
            coarse_expr += self._restrict(self._coarsen_idx(expr, offsets))
        # loop over all terms in the expanded expression and shift the indices
        # of the coarse fields
        r = 0
        for arg in expand(coarse_expr).args:
            shift = self._min_phi_idx(arg)
            if any(shift == np.infty):
                r += arg
            else:
                r += self._shift_index(arg, shift)
        return r

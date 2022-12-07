from itertools import product
import numpy as np
from sympy import Indexed, Symbol, Add, Mul, Pow, expand


class Action:
    """Class representing an action

    An action is represented by an expression which is the sum
    of terms of the form

    F(I)*phi_{I+\alpha}*phi_{I+\beta}

    for d-dimensional summation multiindices I = (i_0,...,i_{d-1}) and
    shift vectors \alpha, \beta. F(I) is some other expression not
    involving phi, which can depend on the summation indices.
    """

    def __init__(self, expr, phi, sidxs):
        """Create new instance

        :arg expr: expression describing the action
        :arg phi: field variable phi
        :arg sidxs: summation multiindices I = (i_0,...,i_{d-1})
        """
        self.phi = phi
        self.sidxs = sidxs
        self.dim = len(sidxs)
        assert self._is_valid_expr(expr), "Not a valid expression " + str(expr)
        self.expr = expand(expr)

    def _get_factors(self, expr):
        """Get a list of factors in a given product expression, separated
        into terms containing phi and not containing phi

        Returns a list form

        [ [list of terms containing phi],
          [list of terms not containing phi] ]

        :arg expr:
        """
        factors = self._get_factors_rec(expr)
        factors_phi = list(
            filter(
                lambda factor: factor.base == self.phi,
                filter(lambda factor: type(factor) == Indexed, factors),
            )
        )
        factors_other = [factor for factor in factors if factor not in factors_phi]
        return factors_phi, factors_other

    def _get_factors_rec(self, expr):
        """Recursively extract a list of factors for a given expression

        Returns None if the expression is not a product

        :arg expr: expression to factor
        """
        if expr.func is Mul:
            result = []
            for term in expr.args:
                result += self._get_factors_rec(term)
            return result
        elif expr.func is Pow:
            return expr.args[1] * self._get_factors_rec(expr.args[0])
        else:
            return [expr]

    def _is_valid_term(self, expr):
        """Verify that the expression is of the form

        F(I)*phi_{I+\alpha}*phi_{I+\beta}

        :arg expr: expression to check
        """
        factors = self._get_factors(expr)
        if factors is None:
            return False
        else:
            factors_phi, _ = self._get_factors(expr)
            # Check that the term contains exactly two terms of the
            # form phi_{I+\alpha}
            valid = len(factors_phi) == 2
            # Check that at least one of the terms is phi_{I} (i.e. there is
            # no shift \alpha)
            shift_indices = [
                [idx.subs([(ell, 0) for ell in self.sidxs]) for idx in term.indices]
                for term in factors_phi
            ]
            valid = valid and (self.dim * [0] in shift_indices)
            return valid

    def _is_valid_expr(self, expr):
        """Verify that the entire expression represents a valid action, i.e.
        it can be written as a sum of valid terms

        :arg expr: expression to check
        """
        expanded_expr = expand(expr)
        if self._is_valid_term(expanded_expr):
            return True
        elif expanded_expr.func is Add:
            return all(map(self._is_valid_term, expanded_expr.args))
        else:
            return False

    def _get_shift_indices(self, expr):
        """Extract the shift indices \alpha from a term
        of the form F(I)*phi_{I}*phi_{I+\alpha}.

        Returns list \alpha

        :arg expr: expression to extract the shift indices from
        """
        phi_terms, _ = self._get_factors(expr)
        shift_indices = [
            [idx.subs([(ell, 0) for ell in self.sidxs]) for idx in term.indices]
            for term in phi_terms
        ]
        if shift_indices[0] == self.dim * [0]:
            return np.asarray(shift_indices[1])
        else:
            return np.asarray(shift_indices[0])

    def _shift_indices(self, x, shift_idxs):
        """Shift the indices of the Indexed variable x by shift_idxs.
        Leave all other expressions unchanged"""
        if type(x) is Indexed:
            idxs = [
                idx.subs(
                    [(ell, ell + shift_idxs[j]) for j, ell in enumerate(self.sidxs)]
                )
                for idx in x.indices
            ]
            return x.base[idxs]
        else:
            return x

    def _mul_expr(self, expr_a, expr_b):
        """Multiply two expressions

        :arg expr_a: first expresion to multiply
        :arg expr_b: second expresion to multiply
        """
        if expr_a.func is Add:
            # Recursive expansion of summands in the first expression
            return sum((self._mul_expr(arg, expr_b) for arg in expr_a.args))
        else:
            if expr_b.func is Add:
                # Recursive expansion of summands in the second expression
                return sum((self._mul_expr(expr_a, arg) for arg in expr_b.args))
            else:
                # Multiply two individual terms
                idxs_a = self._get_shift_indices(expr_a)
                idxs_b = self._get_shift_indices(expr_b)
                # Extract all factors that do not contain phi from factors_a
                _, coeffs_a = self._get_factors(expr_a)
                _, coeffs_b = self._get_factors(expr_b)
                factor = 1
                for coeff in coeffs_a:
                    factor *= coeff
                for coeff in coeffs_b:
                    factor *= self._shift_indices(coeff, idxs_a)
                factor *= self.phi[self.sidxs]
                factor *= self._shift_indices(self.phi[self.sidxs], idxs_a + idxs_b)
                return factor

    def __mul__(self, other):
        """Multiply one action by a different action or an integer

        :arg other: other action to multiply by
        """
        if type(other) is type(self):
            return Action(self._mul_expr(self.expr, other.expr), self.phi, self.sidxs)
        else:
            return Action(other * self.expr, self.phi, self.sidxs)

    def __pow__(self, n):
        """Raise action to n-th integer power

        :arg n: power, must be an integer
        """
        assert type(n) is int
        assert n > 0
        expr = self.expr
        for _ in range(n - 1):
            expr = self._mul_expr(self.expr, expr)
        return Action(expr, self.phi, self.sidxs)


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

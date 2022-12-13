from sympy import Indexed, Add, Mul, Pow, Idx, ccode


class IndexFlattener2D:
    def __init__(self, Mx_lat="Mx_lat", My_lat="My_lat"):
        self.Mx_lat = Mx_lat
        self.My_lat = My_lat

    def _build_hash(self, expr):
        cls = expr.func
        if cls is Add or cls is Mul or cls is Pow:
            # for Add, Mul or Pow return minimum over all subexpressions
            return cls(*[self._build_hash(arg) for arg in expr.args])
        elif cls is Indexed:
            idxs = expr.indices
            idx_hash = "k{HASH" + str(hash(idxs)) + "HASH}"
            self.replacement_dict[
                idx_hash
            ] = f"(({idxs[0]})+{self.Mx_lat:s}) % {self.Mx_lat:s} + {self.Mx_lat:s}*((({idxs[1]})+{self.My_lat:s}) % {self.My_lat:s})"
            k = Idx(idx_hash, 1)
            return expr.base[k]
        else:
            return expr

    def apply(self, expr):
        self.replacement_dict = dict()
        s = ccode(self._build_hash(expr))
        for key, value in self.replacement_dict.items():
            s = s.replace(key, value)
        return s

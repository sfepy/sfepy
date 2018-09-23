# TODO Implement equation in manner similiar to sfepy - it will take of calLing terms evaulation and
# then assembling the matrix, or maybe use the equation from sfepy? hm?

import numpy as nm
from dg_terms import Term


class Equation:
    def __init__(self, terms):
        self.terms = terms
        if len(terms) > 0:
            self.mesh = terms[0].mesh


    def assemble(self):

        # TODO how to get number of cells from mesh?
        A = nm.zeros((len(self.mesh.coors), len(self.mesh.coors)), dtype=nm.float64)

        for term in self.terms:
            # TODO pass arguments to terms?
            val, iels = term.evaluate()
            #from sfepy.base.base import debug; debug()
            A[iels] = A[iels] + val

        return A


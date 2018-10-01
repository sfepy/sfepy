# TODO Implement equation in manner similiar to sfepy - it will take of calLing terms evaulation and
# then assembling the matrix, or maybe use the equation from sfepy? hm?

import numpy as nm
from dg_terms import Term


class Equation:

    def __init__(self, terms):
        self.terms = terms
        if len(terms) > 0:
            self.mesh = terms[0].mesh

    def evaluate(self, mode="weak", dw_mode="vector", asm_obj=None, diff_var=None):

        for term in self.terms:
            val, iels = term.evaluate(diff_var=diff_var)
            term.assemble_to(asm_obj, val, iels, mode=dw_mode)

        return asm_obj





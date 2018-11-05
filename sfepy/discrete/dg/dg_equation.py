import numpy as nm
from dg_terms import DGTerm


class Equation:

    def __init__(self, terms):
        self.terms = terms
        if len(terms) > 0:
            self.mesh = terms[0].mesh

    def evaluate(self, mode="weak", dw_mode="vector", asm_obj=None, diff_var=None, **kwargs):

        for term in self.terms:
            val, iels = term.evaluate(mode=mode, diff_var=diff_var, standalone=False, **kwargs)
            term.assemble_to(asm_obj, val, iels, mode=dw_mode)

        return asm_obj





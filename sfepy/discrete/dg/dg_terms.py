import numpy as nm


class Term:

    def __init__(self, mesh):
        self.mesh = mesh

    def evaluate(self, a):
        raise NotImplemented

    def assemble_to(self, asm_obj, val, iels, mode="vector"):
        if (asm_obj is not None) and (iels is not None):
            asm_obj[iels] = asm_obj[iels] + val


class AdvIntTerm(Term):

    def __init__(self, mesh, a):
        Term.__init__(self, mesh)
        self.a = a
        self.vvar = "v"
        self.diff_var = "u"

    def evaluate(self, mode="weak", diff_var="u"):
        if diff_var == self.diff_var:
            # so far only for approx of order 0
            val = (self.mesh.coors[1:] - self.mesh.coors[:-1]).T  # integral over element with constant test
            # function is just volume of the element
            iels = (nm.arange(len(self.mesh.coors) - 1), nm.arange(len(self.mesh.coors) - 1))
            # values go on to the diagonal, in sfepy this is assured
            # by mesh connectivity induced by basis
        else:
            val = None
            iels = None

        return val, iels


class AdvFluxTerm(Term):

    def __init__(self, mesh, a):
        Term.__init__(self, mesh)
        self.a = a
        self.vvar = None
        self.diff_var = None

    def evaluate(self, mode="weak", diff_var=None, *args):
        if diff_var == self.diff_var:
            # so far only for approx of order 0
            val = nm.ones((len(self.mesh.coors)-1, 1)) * self.a + \
                  nm.ones((len(self.mesh.coors) - 1, 1)) * self.a
            # TODO values of flux terms are functions of solution on previous time step,
            # how to pas these values to the term?
            iels = nm.arange(len(self.mesh.coors) - 1)  # just fill the vector
        else:
            val = None
            iels = None


        return val, iels

import numpy as nm
from sfepy.terms.terms import Term

class DGTerm:

    def __init__(self, mesh):
        self.mesh = mesh

    def evaluate(self, mode="weak", diff_var="u",
                 standalone=True, ret_status=False, **kwargs):
        raise NotImplemented

    def assemble_to(self, asm_obj, val, iels, mode="vector"):
        if (asm_obj is not None) and (iels is not None):
            if mode == "vector":
                if (len(iels) == 2) and (nm.shape(val)[0] == len(iels[0])):
                    for ii in iels[0]:
                        asm_obj[ii][iels[1]] = (asm_obj[ii][iels[1]].T + val[ii]).T
                else:
                    asm_obj[iels] = asm_obj[iels] + val

            elif mode == "matrix":
                if (len(iels) == 3) and (nm.shape(val)[0] == len(iels[0])):
                    for ii in iels[0]:
                        asm_obj[ii][iels[1], iels[2]] = asm_obj[ii][iels[1], iels[2]] + val[ii]
                else:
                    asm_obj[iels] = asm_obj[iels] + val
            else:
                raise ValueError("Unknown assebmly mode '%s'" % mode)


class AdvIntDGTerm(Term):
    # TODO try inheritigng directly from Term?
    def __init__(self, mesh):
        DGTerm.__init__(self, mesh)
        self.vvar = "v"
        self.diff_var = "u"

    def get_fargs(self, *args, **kwargs):

        val = nm.vstack(((self.mesh.coors[1:] - self.mesh.coors[:-1]).T,
                         (self.mesh.coors[1:] - self.mesh.coors[:-1]).T/3))
        # integral over element with constant test
        # function is just volume of the element
        iels = ([0, 1], nm.arange(len(self.mesh.coors) - 1), nm.arange(len(self.mesh.coors) - 1))
        # values go on to the diagonal, in sfepy this is assured
        # by mesh connectivity induced by basis
        fargs = (val, iels)
        return fargs


    def function(self, out, *fargs):

        status = None
        return status

class AdvFluxDGTerm(Term):

    def __init__(self, mesh, a):
        DGTerm.__init__(self, mesh)
        self.a = a
        self.vvar = None
        self.diff_var = None

    def get_fargs(self, state, mode="weak", diff_var=None,
                 standalone=True, ret_status=False, **kwargs):

        u = self.get(state, 'val', step=-1)
        if diff_var == self.diff_var:
            # exact integral
            intg = self.a * (u[0, 1:-1] * (self.mesh.coors[1:] - self.mesh.coors[:-1]) +
                             1/2*u[1, 1:-1]**2 * (self.mesh.coors[1:] - self.mesh.coors[:-1])).T

            fp = self.a * u[0, 2:].T if self.a > 0 else self.a * u[0, 1:-1].T
            fl = self.a * u[0, 1:-1].T if self.a > 0 else self.a * u[0, :-2].T

            val = nm.vstack((fl - fp, - fl - fp + intg))
            # TODO values of flux terms are functions of solution on previous time step,
            # how to pas these values to the term?

            # placement is simple, bud getting the values requires looping over neighbours
            iels = ([0, 1], nm.arange(len(self.mesh.coors) - 1))  # just fill the vector
        else:
            val = None
            iels = None
        fargs = (val, iels)
        return fargs

    def function(selfself, out, *fargs):


        status = None
        return status

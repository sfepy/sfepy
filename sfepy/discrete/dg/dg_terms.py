import numpy as nm


class Term:

    def __init__(self, mesh):
        self.mesh = mesh

    def evaluate(self, a):
        raise NotImplemented


class AdvIntTerm(Term):

    def __init__(self, mesh, a):
        Term.__init__(self, mesh)
        self.a = a


    def evaluate(self):
        # so far only for approx of order 0
        val = (self.mesh.coors[1:] - self.mesh.coors[:-1]).T  # integral over element with constant test
        # function is just volume of the element
        iels = (nm.arange(len(self.mesh.coors) - 1), nm.arange(len(self.mesh.coors) - 1))
        # values go on to the diagonal,  in sfepy this is assured
        # by mesh connectivity induced by basis

        return val, iels


class AdvFluxTerm(Term):

    def __init__(self, mesh, a):
        Term.__init__(self, mesh)
        self.a = a

    def evaluate(self):
        # so far only for approx of order 0
        val = nm.ones((1, 2*(len(self.mesh.coors)-1))) * self.a
        # TODO values of flux terms are functions of solution, how does this go into the matrix?
        iels = (nm.hstack((nm.arange(len(self.mesh.coors) - 1),
                           nm.arange(len(self.mesh.coors) - 1))),
                nm.hstack((nm.arange(len(self.mesh.coors) - 1),
                           nm.arange(len(self.mesh.coors) - 1) + 1)))
        # they should go onto diagonal and one consecutive column

        return val, iels

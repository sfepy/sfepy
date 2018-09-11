import numpy as nm
from sfepy.discrete.fem import Mesh


class Term:

    def __init__(self, mesh):
        self.mesh = mesh


    def evaluate(self, a):
        raise NotImplemented

class AdvIntTerm(Term):

    def evaluate(self, a):


        val = nm.zeros((2, 2))
        iels = nm.ones((2, 2))

        return val, iels


class AdvFluxTerm(Term):


    def evaluate(self):
        val = nm.zeros((2, 2))
        iels = nm.ones((2, 2))

        return val, iels

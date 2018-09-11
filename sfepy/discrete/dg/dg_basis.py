import numpy as nm
# TODO Create basis function for reference elemnt, where does mapping come from? maybe try using sfepy one?
class DGBasis:

    def __init__(self, degree):
        self.degree = degree
        self, degree
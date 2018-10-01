

import numpy as nm
from numpy import dot
import matplotlib.pyplot as plt

# TODO Probably implement Runge-Kutta - will it be possible to separate time stepping solver from
# discretization - it should be ;-)


class TSSolver:

    def __init__(self, eq, ic):
        self.equation = eq
        self.mesh = eq.mesh
        self.initial_cond = ic

    def solve(self, t0, tend, tsteps=10):
        print("Running testing solver: it only calculates matrix A and vector b")
        # TODO how to get number of cells from mesh?
        A = nm.zeros((len(self.mesh.coors), len(self.mesh.coors)), dtype=nm.float64)
        b = nm.zeros((len(self.mesh.coors), 1), dtype=nm.float64)
        A = self.equation.evaluate(dw_mode="matrix", asm_obj=A, diff_var="v")
        self.equation.evaluate(dw_mode="vector", asm_obj=b, diff_var=None)
        plt.imshow(A)
        plt.show()

        # U = dot(nm.linalg.inv(A), b)

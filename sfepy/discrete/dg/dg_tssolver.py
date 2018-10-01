

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
        print("Running testing solver: it does not solve anything, only tests shape and types of data!")
        # TODO how to get number of cells from mesh?
        A = nm.zeros((len(self.mesh.coors)-1, len(self.mesh.coors)-1), dtype=nm.float64)
        b = nm.zeros((len(self.mesh.coors)-1, 1), dtype=nm.float64)
        u = nm.ones((len(self.mesh.coors), 1), dtype=nm.float64)
        u[5] = 2
        self.equation.evaluate(dw_mode="matrix", asm_obj=A, diff_var="u")
        self.equation.evaluate(dw_mode="vector", asm_obj=b, diff_var=None, u=u)

        print(A)
        print(b)

        plt.imshow(A)
        plt.plot(b)
        plt.show()

        U = dot(nm.linalg.inv(A), b)

        print(U)
        plt.plot(U)
        plt.show()

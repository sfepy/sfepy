

import numpy as nm
from numpy import dot
import matplotlib.pyplot as plt

# TODO Probably implement Runge-Kutta - will it be possible to separate time stepping solver from
# discretization - it should be ;-)


class TSSolver:

    def __init__(self, eq, ic, bc):
        self.equation = eq
        self.mesh = eq.mesh
        self.initial_cond = ic
        self.boundary_cond = bc

    def solve(self, t0, tend, tsteps=10):
        print("Running testing solver: it does not solve anything, only tests shapes and types of data!")
        # TODO how to get number of cells from mesh?
        A = nm.zeros((2, len(self.mesh.coors)-1, len(self.mesh.coors)-1), dtype=nm.float64)
        b = nm.zeros((2, len(self.mesh.coors)-1, 1), dtype=nm.float64)
        u = nm.ones((2, len(self.mesh.coors) + 1 , 1), dtype=nm.float64)

        # ic
        u[0, 1:-1] = self.initial_cond
        u[1, 1:-1] = self.initial_cond

        # bc
        u[0, 0] = self.boundary_cond["left"]
        u[0, -1] = self.boundary_cond["right"]
        u[1, 0] = self.boundary_cond["left"]
        u[1, -1] = self.boundary_cond["right"]

        self.equation.evaluate(dw_mode="matrix", asm_obj=A, diff_var="u")
        self.equation.evaluate(dw_mode="vector", asm_obj=b, diff_var=None, u=u)

        # print(A)
        # print(b)

        plt.figure("A[0]")
        plt.imshow(A[0])
        plt.figure("A[1]")
        plt.imshow(A[1])

        plt.figure("b")
        plt.plot(b[0], label="b0")
        plt.plot(b[1], label="b1")
        plt.legend()

        u[0, 1:-1] = dot(nm.linalg.inv(A[0]), b[0])
        u[1, 1:-1] = dot(nm.linalg.inv(A[1]), b[1])

        # print(u)
        plt.figure("u")
        plt.plot(u[0], label="u0")
        plt.plot(u[1], label="u1")
        plt.legend()
        plt.show()

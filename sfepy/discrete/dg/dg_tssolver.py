

import numpy as nm
from numpy import dot
import matplotlib.pyplot as plt


class TSSolver:

    def __init__(self, eq, ic, bc):
        self.equation = eq
        self.mesh = eq.mesh
        # TODO setup initial condition - integrate over elements and set coefs
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

class RK3Solver(TSSolver):

    def solve(self, t0, tend, tsteps=10):
        # TODO implement three step Ruge-Kutta
        dt = (tend - t0) / tsteps

        A = nm.zeros((2, len(self.mesh.coors) - 1, len(self.mesh.coors) - 1), dtype=nm.float64)
        b = nm.zeros((2, len(self.mesh.coors) - 1, 1), dtype=nm.float64)
        u = nm.ones((2, len(self.mesh.coors) + 1, tsteps), dtype=nm.float64)
        u1 = nm.ones((2, len(self.mesh.coors) + 1), dtype=nm.float64)
        u2 = nm.ones((2, len(self.mesh.coors) + 1), dtype=nm.float64)
        u3 = nm.ones((2, len(self.mesh.coors) + 1), dtype=nm.float64)

        # TODO check shapes and simplify bc treatment
        # bc
        u[0, 0, 0] = self.boundary_cond["left"]
        u[0, -1, 0] = self.boundary_cond["right"]
        u[1, 0, 0] = self.boundary_cond["left"]
        u[1, -1, 0] = self.boundary_cond["right"]

        # ic
        u[0, 1:-1, 0] = self.initial_cond[0, :]
        u[1, 1:-1, 0] = self.initial_cond[1, :]

        for it in range(1, tsteps):
            # 1st stage
            # bcs
            u1[0, 0] = self.boundary_cond["left"]
            u1[0, -1] = self.boundary_cond["right"]
            u1[1, 0] = self.boundary_cond["left"]
            u1[1, -1] = self.boundary_cond["right"]

            # get RHS
            self.equation.evaluate(dw_mode="matrix", asm_obj=A, diff_var="u")
            self.equation.evaluate(dw_mode="vector", asm_obj=b, diff_var=None, u=u[:, :, it-1])

            # get update u1
            u1[0, 1:-1] = u[0, 1:-1, it-1] + dt * dot(nm.linalg.inv(A[0]), b[0])
            u1[1, 1:-1] = u[1, 1:-1, it-1] + dt * dot(nm.linalg.inv(A[1]), b[1])

            # 2nd stage
            # bcs
            u2[0, 0] = self.boundary_cond["left"]
            u2[0, -1] = self.boundary_cond["right"]
            u2[1, 0] = self.boundary_cond["left"]
            u2[1, -1] = self.boundary_cond["right"]

            # get RHS
            self.equation.evaluate(dw_mode="matrix", asm_obj=A, diff_var="u")
            self.equation.evaluate(dw_mode="vector", asm_obj=b, diff_var=None, u=u1[:, :])

            # get update u2
            u2[0, 1:-1] = 3 * u[0, 1:-1, it - 1]/4  \
                            + u1[0, 1:-1] / 4 \
                            + dt * dot(nm.linalg.inv(A[0]), b[0]) / 4
            u2[1, 1:-1] = 3 * u[1, 1:-1, it - 1]/4  \
                            + u1[1, 1:-1] / 4 \
                            + dt * dot(nm.linalg.inv(A[1]), b[1]) / 4

            # 3rd stage
            # bcs
            u3[0, 0] = self.boundary_cond["left"]
            u3[0, -1] = self.boundary_cond["right"]
            u3[1, 0] = self.boundary_cond["left"]
            u3[1, -1] = self.boundary_cond["right"]

            # get RHS
            self.equation.evaluate(dw_mode="matrix", asm_obj=A, diff_var="u")
            self.equation.evaluate(dw_mode="vector", asm_obj=b, diff_var=None, u=u2[:, :])

            # get update u3
            u3[0, 1:-1] = u[0, 1:-1, it - 1] / 3 \
                          + u2[0, 1:-1] / 3 \
                          + dt * dot(nm.linalg.inv(A[0]), b[0]) / 3
            u3[1, 1:-1] = u[1, 1:-1, it - 1] / 3 \
                          + u2[1, 1:-1] / 3 \
                          + dt * dot(nm.linalg.inv(A[1]), b[1]) / 3

            u[:, :, it-1] = u3


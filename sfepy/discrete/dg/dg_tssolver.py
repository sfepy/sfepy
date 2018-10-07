

import numpy as nm
from numpy import dot
import matplotlib.pyplot as plt

from numpy import newaxis as nax

class TSSolver:

    def __init__(self, eq, ic, bc, basis):
        self.equation = eq
        self.mesh = eq.mesh
        self.basis = basis
        self.initial_cond = self.sampleIC(self.mesh, ic, self.intGauss2, self.basis)
        self.boundary_cond = bc

    def sampleIC(self, mesh, ic, quad, basis):

        sic = nm.zeros((2, len(self.mesh.coors)-1, 1), dtype=nm.float64)
        sic[0, :] = quad(mesh, ic)/2
        sic[1, :] = 2*quad(mesh, ic)/3

        return sic

    @staticmethod
    def intGauss2(mesh, f):
        # TODO check transformation to the reference element
        x_1 = mesh.coors[:-1] + (mesh.coors[1:] - mesh.coors[:-1]) * (-nm.sqrt(1 / 2) + 1) / 2
        x_2 = mesh.coors[:-1] + (mesh.coors[1:] - mesh.coors[:-1]) * (nm.sqrt(1 / 3) + 1) / 2

        w = (mesh.coors[1:] - mesh.coors[:-1]) / 2

        return f(x_1) +  f(x_2)

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
        dt = (tend - t0) / tsteps

        A = nm.zeros((2, len(self.mesh.coors) - 1, len(self.mesh.coors) - 1), dtype=nm.float64)
        b = nm.zeros((2, len(self.mesh.coors) - 1, 1), dtype=nm.float64)
        u  = nm.ones((2, len(self.mesh.coors) + 1, tsteps, 1), dtype=nm.float64)
        u1 = nm.ones((2, len(self.mesh.coors) + 1, 1), dtype=nm.float64)
        u2 = nm.ones((2, len(self.mesh.coors) + 1, 1), dtype=nm.float64)
        u3 = nm.ones((2, len(self.mesh.coors) + 1, 1), dtype=nm.float64)

        # bc
        u[:, 0, 0] = self.boundary_cond["left"]
        u[:, -1, 0] = self.boundary_cond["right"]

        # ic
        u[:, 1:-1, 0] = self.initial_cond

        for it in range(1, tsteps):
            # ----1st stage----
            # bcs
            u1[:, 0] = self.boundary_cond["left"]
            u1[:, -1] = self.boundary_cond["right"]

            # get RHS
            self.equation.evaluate(dw_mode="matrix", asm_obj=A, diff_var="u")
            self.equation.evaluate(dw_mode="vector", asm_obj=b, diff_var=None, u=u[:, :, it-1])

            # get update u1
            u1[0, 1:-1] = u[0, 1:-1, it-1] + dt * dot(nm.linalg.inv(A[0]), b[0])
            u1[1, 1:-1] = u[1, 1:-1, it-1] + dt * dot(nm.linalg.inv(A[1]), b[1])

            # ----2nd stage----
            # bcs
            u2[:, 0] = self.boundary_cond["left"]
            u2[:, -1] = self.boundary_cond["right"]

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

            # ----3rd stage-----
            # bcs
            u3[:, 0] = self.boundary_cond["left"]
            u3[:, -1] = self.boundary_cond["right"]

            # get RHS
            self.equation.evaluate(dw_mode="matrix", asm_obj=A, diff_var="u")
            self.equation.evaluate(dw_mode="vector", asm_obj=b, diff_var=None, u=u2[:, :])

            # get update u3
            u3[0, 1:-1] = u[0, 1:-1, it - 1] / 3 \
                          + 2*u2[0, 1:-1] / 3 \
                          + 2*dt * dot(nm.linalg.inv(A[0]), b[0]) / 3
            u3[1, 1:-1] = u[1, 1:-1, it - 1] / 3 \
                          + 2*u2[1, 1:-1] / 3 \
                          + 2*dt * dot(nm.linalg.inv(A[1]), b[1]) / 3

            u[:, :, it] = u3

        return u, dt


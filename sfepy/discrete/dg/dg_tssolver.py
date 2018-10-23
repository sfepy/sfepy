

import numpy as nm
from numpy import dot
import matplotlib.pyplot as plt

from numpy import newaxis as nax


class TSSolver:
    # TODO refactor Solver to Problem class

    def __init__(self, eq, ic, bc, basis):
        self.equation = eq
        self.mesh = eq.mesh
        self.basis = basis
        self.initial_cond = self.sampleIC(self.mesh, ic, self.intGauss2, self.basis)
        self.boundary_cond = bc

    def sampleIC(self, mesh, ic, quad, basis):
        sic = nm.zeros((2, self.mesh.n_el, 1), dtype=nm.float64)

        # FIXME initialization still is not correct!
        # TODO check transformation to the reference element
        c = (mesh.coors[1:] + mesh.coors[:-1])/2  # center
        s = (mesh.coors[1:] - mesh.coors[:-1])/2  # scale
        sic[0, :] = quad(lambda t: ic(c + t*s))/2
        sic[1, :] = 3*quad(lambda t: t*ic(c + t*s))/2
        return sic

    @staticmethod
    def intGauss2(f):

        x_1 = -1/nm.sqrt(1./3)
        x_2 = 1/nm.sqrt(1./3)

        return f(x_1) + f(x_2)

    @staticmethod
    def intGauss3(f):
        x_0 = 0
        x_1 = - nm.sqrt(3./5.)
        x_2 = nm.sqrt(3./5.)

        w_0 = 8./9.
        w_1 = 5. / 9.

        return w_0 * f(x_0) + w_1 * f(x_1) + w_1 * f(x_2)

    def solve(self, t0, tend, tsteps=10):
        print("Running testing solver: it does not solve anything, only tests shapes and types of data!")
        A = nm.zeros((2, self.mesh.n_el, self.mesh.n_el), dtype=nm.float64)
        b = nm.zeros((2, self.mesh.n_el-1, 1), dtype=nm.float64)
        u = nm.zeros((2, self.mesh.n_nod + 1 , 1), dtype=nm.float64)

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
        dt = float(tend - t0) / tsteps
        dx = nm.max(self.mesh.coors[1:] - self.mesh.coors[:-1])
        dtdx = dt/dx
        maxa = abs(self.equation.terms[1].a)

        print("Space divided into {0} cells, {1} steps, step size is {2}".format(self.mesh.n_el, len(self.mesh.coors), dx))
        print("Time divided into {0} nodes, {1} steps, step size is {2}".format(tsteps - 1, tsteps, dt))
        print("Courant number c = max(abs(u)) * dt/dx = {0}".format(maxa * dtdx))

        A  = nm.zeros((2, self.mesh.n_el, self.mesh.n_el), dtype=nm.float64)
        b  = nm.zeros((2, self.mesh.n_el, 1), dtype=nm.float64)
        u  = nm.zeros((2, self.mesh.n_el + 2, tsteps, 1), dtype=nm.float64)
        u1 = nm.zeros((2, self.mesh.n_el + 2, 1), dtype=nm.float64)
        u2 = nm.zeros((2, self.mesh.n_el + 2, 1), dtype=nm.float64)

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
            A[:] = 0
            b[:] = 0
            self.equation.evaluate(dw_mode="matrix", asm_obj=A, diff_var="u")
            self.equation.evaluate(dw_mode="vector", asm_obj=b, diff_var=None, u=u[:, :, it-1])

            # get update u1
            u1[0, 1:-1] = u[0, 1:-1, it-1] + dt * b[0] / nm.diag(A[0])[:, nax]
            u1[1, 1:-1] = u[1, 1:-1, it-1] + dt * b[1] / nm.diag(A[1])[:, nax]

            # ----2nd stage----
            # bcs
            u2[:, 0] = self.boundary_cond["left"]
            u2[:, -1] = self.boundary_cond["right"]

            # get RHS
            A[:] = 0
            b[:] = 0
            self.equation.evaluate(dw_mode="matrix", asm_obj=A, diff_var="u")
            self.equation.evaluate(dw_mode="vector", asm_obj=b, diff_var=None, u=u1[:, :])

            # get update u2
            u2[0, 1:-1] = (3 * u[0, 1:-1, it - 1] + u1[0, 1:-1]
                           + dt * b[0] / nm.diag(A[0])[:, nax]) / 4
            u2[1, 1:-1] = (3 * u[1, 1:-1, it - 1] + u1[1, 1:-1]
                           + dt * b[1] / nm.diag(A[1])[:, nax]) / 4

            # ----3rd stage-----
            # get RHS
            A[:] = 0
            b[:] = 0
            self.equation.evaluate(dw_mode="matrix", asm_obj=A, diff_var="u")
            self.equation.evaluate(dw_mode="vector", asm_obj=b, diff_var=None, u=u2[:, :])

            # get update u3
            u[0, 1:-1, it] = (u[0, 1:-1, it - 1] + 2 * u2[0, 1:-1]
                              + 2*dt * b[0] / nm.diag(A[0])[:, nax]) / 3
            u[1, 1:-1, it] = (u[1, 1:-1, it - 1] + 2 * u2[1, 1:-1]
                              + 2*dt * b[1] / nm.diag(A[1])[:, nax]) / 3

        return u, dt


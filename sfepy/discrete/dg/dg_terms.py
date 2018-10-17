import numpy as nm


class DGTerm:

    def __init__(self, mesh):
        self.mesh = mesh

    def evaluate(self, mode="weak", diff_var="u",
                 standalone=True, ret_status=False, **kwargs):
        raise NotImplemented

    @staticmethod
    def assemble_to(asm_obj, val, iels, mode="vector"):
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
                raise ValueError("Unknown assembly mode '%s'" % mode)


class AdvIntDGTerm(DGTerm):

    def __init__(self, mesh):
        DGTerm.__init__(self, mesh)
        self.vvar = "v"
        self.diff_var = "u"

    def evaluate(self, mode="weak", diff_var="u",
                 standalone=True, ret_status=False, **kwargs):
        if diff_var == self.diff_var:

            val = nm.vstack(((self.mesh.coors[1:] - self.mesh.coors[:-1]).T,
                             (self.mesh.coors[1:] - self.mesh.coors[:-1]).T/3))

            iels = ([0, 1], nm.arange(self.mesh.n_el),
                            nm.arange(self.mesh.n_el))
            # values go on to the diagonal, in sfepy this is assured
            # by mesh connectivity induced by basis
        else:
            val = None
            iels = None

        return val, iels


class AdvFluxDGTerm(DGTerm):
    """
    So far Lax-Friedrichs flux
    """

    def __init__(self, mesh, a):
        DGTerm.__init__(self, mesh)
        self.a = a
        self.vvar = None
        self.diff_var = None

    def evaluate(self, mode="weak", diff_var=None,
                 standalone=True, ret_status=False, **kwargs):

        u = kwargs.pop('u', None)
        if diff_var == self.diff_var:

            # for Legendre basis integral of higher order
            # functions of the basis is zero,
            # hence we calculate integral
            #
            # int_{j-1/2}^{j+1/2} f(u)dx
            #
            # only from the zero order function, over [-1, 1] - hence the 2
            intg = self.a * u[0, 1:-1].T * 2

            #  the Lax-Friedrichs flux is
            #       F(a, b) = 1/2(f(a) + f(b)) + max(f'(w)) / 2 * (a - b)
            # in our case a and b are values to the left and right of the element boundary
            # for Legendre basis these are:
            # u_left = U_0 + U_1 + U_2 + ...
            # u_right = U_0 - U_1 + U_2 + ... = sum_0^{order} (-1)^p * U_p

            # left flux is calculated in j_-1/2  where U(j-1) and U(j) meet
            # right flux is calculated in j_+1/2 where U(j) and U(j+1) meet
            fl = self.a * (u[0, :-2] + u[1, :-2] +
                          (u[0, 1:-1] - u[1, 1:-1])).T / 2 + \
                 nm.abs(self.a) * (u[0, :-2] + u[1, :-2] -
                                  (u[0, 1:-1] - u[1, 1:-1])).T / 2

            fp = self.a * (u[0, 1:-1] + u[1, 1:-1] +
                          (u[0, 2:] - u[1, 2:])).T / 2 + \
                 nm.abs(self.a) * (u[0, 1:-1] + u[1, 1:-1] -
                                  (u[0, 2:] - u[1, 2:])).T / 2

            val = nm.vstack((fl - fp, - fl - fp + intg))

            # placement is simple, but getting the values requires looping over neighbours
            iels = ([0, 1], nm.arange(len(self.mesh.coors) - 1))  # just fill the vector
        else:
            val = None
            iels = None

        return val, iels

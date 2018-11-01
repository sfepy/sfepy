import numpy as nm
from sfepy.terms.terms import Term

class DGTerm:

    def __init__(self, mesh):
        self.mesh = mesh

    def evaluate(self, mode="weak", diff_var="u",
                 standalone=True, ret_status=False, **kwargs):
        raise NotImplemented

    def assemble_to(self, asm_obj, val, iels, mode="vector"):
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
                raise ValueError("Unknown assebmly mode '%s'" % mode)


class AdvIntDGTerm(Term):
    # TODO try inheritigng directly from Term?
    # TODO Replace this term by sfepy.terms.dw_volume?
    name = "dw_volume"

    def get_fargs(self, *args, **kwargs):

        val = nm.vstack(((self.mesh.coors[1:] - self.mesh.coors[:-1]).T,
                         (self.mesh.coors[1:] - self.mesh.coors[:-1]).T/3))
        # integral over element with constant test
        # function is just volume of the element
        iels = ([0, 1], nm.arange(len(self.mesh.coors) - 1), nm.arange(len(self.mesh.coors) - 1))
        # values go on to the diagonal, in sfepy this is assured
        # by mesh connectivity induced by basis
        fargs = (val, iels)
        return fargs


    def function(self, out, *fargs):

        status = None
        return status

class AdvFluxDGTerm(Term):

    name = "dw_dg_advect__flux"
    modes = ("weak",)
    arg_types = ('material', 'virtual', 'state')
    arg_shapes = {'material': 'a, 1', 'virtual': ('1', 'state'),
                  'state'   : '1'}
    symbolic = {'expression' : 'grad(a*u)',
                'map': {'u': 'state', 'a': 'material'}
    }

    def get_fargs(self, state, mode="weak", diff_var=None,
                 standalone=True, ret_status=False, **kwargs):

        if diff_var == self.diff_var:
            u = self.get(state, 'val', step=-1)
            a = self.get()

        else:
            val = None
            iels = None
        fargs = (u, a)
        return fargs

    def function(self, out, u, a):
        # for Legendre basis integral of higher order
        # functions of the basis is zero,
        # hence we calculate integral
        #
        # int_{j-1/2}^{j+1/2} f(u)dx
        #
        # only from the zero order function, over [-1, 1] - hence the 2
        intg = a * u[0, 1:-1].T * 2
        
        #  the Lax-Friedrichs flux is
        #       F(a, b) = 1/2(f(a) + f(b)) + max(f'(w)) / 2 * (a - b)
        # in our case a and b are values to the left and right of the element boundary
        # for Legendre basis these are:
        # u_left = U_0 + U_1 + U_2 + ...
        # u_right = U_0 - U_1 + U_2 + ... = sum_0^{order} (-1)^p * U_p

        # left flux is calculated in j_-1/2  where U(j-1) and U(j) meet
        # right flux is calculated in j_+1/2 where U(j) and U(j+1) meet

        fl = a * (u[0, :-2] + u[1, :-2] +
                 (u[0, 1:-1] - u[1, 1:-1])).T / 2 + \
             nm.abs(a) * (u[0, :-2] + u[1, :-2] -
                         (u[0, 1:-1] - u[1, 1:-1])).T / 2

        fp = a * (u[0, 1:-1] + u[1, 1:-1] +
                 (u[0, 2:] - u[1, 2:])).T / 2 + \
             nm.abs(a) * (u[0, 1:-1] + u[1, 1:-1] -
                         (u[0, 2:] - u[1, 2:])).T / 2

        val = nm.vstack((fl - fp, - fl - fp + intg))

        iels = ([0, 1], nm.arange(len(self.mesh.coors) - 1))  # just fill the vector

        vals = nm.vstack((fl - fp, - fl - fp + intg))
        out = (vals, iels)
        status = None
        return status

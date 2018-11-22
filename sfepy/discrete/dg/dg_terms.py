import numpy as nm
from sfepy.terms.terms import Term


class AdvVolDGTerm(Term):

    name = "dw_dg_volume"
    modes = ("weak",)
    arg_types = ('virtual', 'state')
    arg_shapes = {'virtual': ('D', 'state'),
                  'state': 1}
    symbolic = {'expression': '1',
                'map': {'u': 'state'}
                }

    def __init__(self, integral, region, u=None, v=None):
        Term.__init__(self, "adv_vol(v, u)", "v, u", integral, region, u=u, v=v)
        self.u = u
        self.v = v
        self.setup()

    def get_fargs(self, test, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        if diff_var is not None:
            doeval = True
        else:
            doeval = False

        return (doeval,)

    def function(self, out, doeval):
        if doeval:
            vols = self.region.domain.cmesh.get_volumes(1)
            # TODO which dimension do we really want?
            # integral over element with constant test
            # function is just volume of the element
            out[:] = 0
            # out[:, 0, 0, 0] = vols
            # out[:, 0, 1, 1] = vols / 3.0
            out[::2, 0, 0, 0] = vols
            out[1::2, 0, 0, 0] = vols / 3.0
            # TODO move to for cycle to add values for higher order approx
        else:
            out[:] = 0.0
        status = None
        return status


class AdvFluxDGTerm(Term):

    def __init__(self, integral, region, u=None, v=None, a=lambda x: 1):
        Term.__init__(self, "adv_lf_flux(a.val, v, u)", "a.val, v, u", integral, region, u=u, v=v, a=a)
        self.u = u
        self.v = v
        self.a = a
        self.setup()

    name = "dw_dg_advect_flux"
    modes = ("weak",)
    arg_types = ('material', 'virtual', 'state')
    arg_shapes = {'material': 'D, 1',
                  'virtual': ('D', 'state'),
                  'state'   : 1}
    symbolic = {'expression' : 'grad(a*u)',
                'map': {'u': 'state', 'a': 'material'}
    }

    def get_fargs(self, a, test, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):

        # varc = self.get_variables(as_list=False)['u']
        u = self.get(state, 'dg', step=-1)

        if diff_var is not None:
            doeval = False
        else:
            doeval = True

        fargs = u, a[:, :, 0, 0], doeval
        return fargs

    def function(self, out, u, a, doeval):
        if not doeval:
            out[:] = 0
            return None

        # for Legendre basis integral of higher order
        # functions of the basis is zero,
        # hence we calculate integral
        #
        # int_{j-1/2}^{j+1/2} f(u)dx
        #
        # only from the zero order function, over [-1, 1] - hence the 2
        intg = a[:, 0] * u[:, 0] * 2

        #  the Lax-Friedrichs flux is
        #       F(a, b) = 1/2(f(a) + f(b)) + max(f'(w)) / 2 * (a - b)
        # in our case a and b are values to the left and right of the element boundary
        # for Legendre basis these are:
        # u_left = U_0 + U_1 + U_2 + ...
        # u_right = U_0 - U_1 + U_2 + ... = sum_0^{order} (-1)^p * U_p

        # left flux is calculated in j_-1/2  where U(j-1) and U(j) meet
        # right flux is calculated in j_+1/2 where U(j) and U(j+1) meet

        # TODO how to treat BCs inside term?
        bc_shape = (1, ) + nm.shape(u)[1:]
        bcl = u[0].reshape(bc_shape)
        bcr = u[-1].reshape(bc_shape)
        ur = nm.concatenate((u[1:], bcr))
        ul = nm.concatenate((bcl, u[:-1]))

        # TODO move to for cycle to add values for higher order approx
        # TODO research general fluxes for higher dimensions

        # fl = a * (u[0, :-2] + u[1, :-2] +
        #           (u[0, 1:-1] - u[1, 1:-1])).T / 2 + \
        #      nm.abs(a) * (u[0, :-2] + u[1, :-2] -
        #                   (u[0, 1:-1] - u[1, 1:-1])).T / 2
        # FIXME u and a will most likely have different shape
        fl = a[:, 0] * (ul[:, 0] + ul[:, 1] +
                        (u[:, 0] - u[:, 1])) / 2 + \
            nm.abs(a[:, 0]) * (ul[:, 0] + ul[:, 1] -
                               (u[:, 0] - u[:, 1])) / 2

        # fp = a * (u[0, 1:-1] + u[1, 1:-1] +
        #           (u[0, 2:] - u[1, 2:])).T / 2 + \
        #      nm.abs(a) * (u[0, 1:-1] + u[1, 1:-1] -
        #                   (u[0, 2:] - u[1, 2:])).T / 2

        fp = a[:, 0] * (u[:, 0] + u[:, 1] +
                        (ur[:, 0] - ur[: , 1])) / 2 + \
            nm.abs(a[:, 0]) * (u[:, 0] + u[:, 1] -
                         (ur[:, 0] - ur[:, 1])) / 2

        # val = nm.vstack((fl - fp, - fl - fp + intg))

        out[:] = 0.0
        # out[:, 0, 0, 0] = (fl - fp)[:, 0, 0]
        # out[:, 0, 1, 0] = (- fl - fp + intg)[:, 0, 0]
        out[::2, 0, 0, 0] = (fl - fp)  # this is how DGField should work
        out[1::2, 0, 0, 0] = (- fl - fp + intg)
        status = None
        return status
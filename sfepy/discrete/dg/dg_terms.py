import numpy as nm
from sfepy.terms.terms import Term


class AdvVolDGTerm(Term):

    name = "dw_dg_volume"
    modes = ("weak",)
    arg_types = ('virtual', 'state')
    # arg_types = ('ts', 'virtual', 'state')
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
            mtx_mode = True
            u = None
        else:
            mtx_mode = False
            ur = state.data[0]  # this uses data provided by solver

            n_cell = self.region.get_n_cells(False)
            u = nm.zeros((n_cell, 2))  # 2 is approx order!
            for i in range(2):
                u[:, i] = ur[n_cell * i: n_cell * (i + 1)]

        # vg, _ = self.get_mapping(state)

        return u, mtx_mode

    def function(self, out, u, mtx_mode):
        vols = self.region.domain.cmesh.get_volumes(1)
        if mtx_mode:
            # TODO which dimension do we really want?

            # integral over element with constant test
            # function is just volume of the element
            out[:] = 0.0
            # out[:, 0, 0, 0] = vols
            # out[:, 0, 1, 1] = vols / 3.0
            out[:, 0, 0, 0] = vols
            out[:, 0, 1, 1] = vols / 3.0
            out[:, 0, 2, 2] = 1
            # TODO move to for cycle to add values for higher order approximations
        else:
            out[:] = 0.0
            out[:, 0, 0, 0] = vols * u[:, 0]
            out[:, 0, 1, 0] = vols/3. * u[:, 1]
            # out[:] = 0
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

        if diff_var is not None:
            # do not eval in matrix mode, we however still need
            # this term to have diff_var in order for it to receive the values
            doeval = False
            return None, None, doeval
        else:
            doeval = True
            ur = state.data[0]  # this uses data provided by solver
            # ur = self.get(state, 'dg', step=-1)  # however this probably too,
            # as they set to variable in equation

            # reshape DOFs vector for convenience in function()
            n_cell = state.field.n_cell
            n_nod = state.field.poly_space.n_nod
            u = nm.zeros((n_cell, n_nod))
            for i in range(n_nod):
                u[:, i] = ur[n_cell * i : n_cell*(i+1)]

            fargs = u, a[:, :, 0, 0], doeval
            return fargs

    def function(self, out, u, a, doeval):
        if not doeval:
            out[:] = 0.0
            return None

        # for Legendre basis integral of higher order
        # functions of the basis is zero,
        # hence we calculate integral
        #
        # int_{j-1/2}^{j+1/2} f(u)dx
        #
        # only from the zero order function, over [-1, 1] - hence the 2
        intg = a[:, 0] * u[:, 0] * 2
        # i.e. intg = a * u0 * reference_el_vol

        #  the Lax-Friedrichs flux is

        #       F(a, b) = 1/2(f(a) + f(b)) + max(f'(w)) / 2 * (a - b)

        # in our case a and b are values in the elements left and right of
        # the respective element boundaries
        # for Legendre basis these are:
        # u_left = U_0 + U_1 + U_2 + ...
        # u_right = U_0 - U_1 + U_2 - U_3 ... = sum_{p=0}^{order} (-1)^p * U_p

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

        fl = a[:, 0] * (ul[:, 0] + ul[:, 1] +
                        (u[:, 0] - u[:, 1])) / 2 + \
            nm.abs(a[:, 0]) * (ul[:, 0] + ul[:, 1] -
                               (u[:, 0] - u[:, 1])) / 2

        fp = a[:, 0] * (u[:, 0] + u[:, 1] +
                        (ur[:, 0] - ur[: , 1])) / 2 + \
            nm.abs(a[:, 0]) * (u[:, 0] + u[:, 1] -
                               (ur[:, 0] - ur[:, 1])) / 2

        out[:] = 0.0
        flux0 = (fl - fp)
        flux1 = (- fl - fp + intg)

        out[:, 0, 0, 0] = -flux0
        out[:, 0, 1, 0] = -flux1

        # compute residual
        # vols = self.region.domain.cmesh.get_volumes(1)
        # out[:, 0, 0, 0] = vols * u[:, 0] - flux0
        # out[:, 0, 1, 0] = vols/3 * u[:, 1] - flux1
        status = None
        return status
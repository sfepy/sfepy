import numpy as nm
from sfepy.terms.terms import Term

from dg_field import get_unraveler, get_raveler



def unravel_sol(state):
    """
    Unpacks solution from flat vector to
    (n_cell, n_el_nod, n_comp)


    :param state:
    """

    u = state.data[0]  # this uses data provided by solver
    # this uses data provided by solver
    # ur = self.get(state, 'dg', step=-1)
    # however this would probably too,
    # as they are set to variable in equation

    # in 2D+ case this will be replaced by get_nbrhd_dofs

    n_cell = state.field.n_cell
    n_el_nod = state.field.poly_space.n_nod
    unravel = get_unraveler(n_el_nod, n_cell)
    return unravel(u)


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
            n_el_nod = state.field.poly_space.n_nod
        else:
            mtx_mode = False
            u = unravel_sol(state)
            n_el_nod = state.field.poly_space.n_nod


        return u, mtx_mode, n_el_nod

    def function(self, out, u, mtx_mode, n_el_nod):
        vols = self.region.domain.cmesh.get_volumes(1)[:, None]
        if mtx_mode:
            # TODO which dimension do we really want?

            # integral over element with constant test
            # function is just volume of the element
            out[:] = 0.0
            # out[:, 0, 0, 0] = vols
            # out[:, 0, 1, 1] = vols / 3.0
            for i in range(n_el_nod):
                out[:, :, i, i] = vols / (2.0 * i + 1.0)
        else:
            out[:] = 0.0
            for i in range(n_el_nod):
                out[:, :, i, 0] = vols / (2.0 * i + 1.0) * u[:, i]
                # TODO does this hold for higher orders? -  seems like it does...
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
            return None, None, None, None, doeval, 0
        else:
            doeval = True

            u = unravel_sol(state) # TODO we unravel twice, refactor
            n_el_nod = state.field.poly_space.n_nod
            nb_dofs, nb_normals = state.field.get_nbrhd_dofs(state.field.region, state)
            # state variable has dt in it!

            fargs = (u, nb_dofs, nb_normals, a[:, :1, 0, 0], doeval, n_el_nod)
            return fargs

    def function(self, out, u, nb_u, nb_n, velo, doeval, n_el_nod):
        if not doeval:
            out[:] = 0.0
            return None

        #  the Lax-Friedrichs flux is

        #       F(a, b) = 1/2(f(a) + f(b)) + max(|f'(w)|) / 2 * (a - b)

        # in our case a and b are values in the elements left and right of
        # the respective element boundaries
        # for Legendre basis these are:
        # u_left = U_0 + U_1 + U_2 + ...
        # u_right = U_0 - U_1 + U_2 - U_3 ... = sum_{p=0}^{order} (-1)^p * U_p

        # left flux is calculated in j_-1/2  where U(j-1) and U(j) meet
        # right flux is calculated in j_+1/2 where U(j) and U(j+1) meet

        # fl:
        # fl = velo[:, 0] * (ul[:, 0] + ul[:, 1] +
        #                    (u[:, 0] - u[:, 1])) / 2 + \
        #      nm.abs(velo[:, 0]) * (ul[:, 0] + ul[:, 1] -
        #                            (u[:, 0] - u[:, 1])) / 2
        a = 0
        b = 0
        sign = 1
        for i in range(n_el_nod):
            a += nb_u[:, 0, i]  # TODO iterate over nuber n_el_facets
            b += sign * u[:, i]
            sign *= -1

        fl = (velo * a + velo * b) / 2 + \
              nm.abs(velo) * (a - b) / 2

        # fl:
        # fp = velo[:, 0] * (u[:, 0] + u[:, 1] +
        #                            (ur[:, 0] - ur[: , 1])) / 2 + \
        #              nm.abs(velo[:, 0]) * (u[:, 0] + u[:, 1] -
        #                                    (ur[:, 0] - ur[:, 1])) / 2
        a = 0
        b = 0
        sign = 1
        for i in range(n_el_nod):
            a += u[:, i]
            b += sign * nb_u[:,1 , i]
            sign *= -1

        fp = (velo * a + velo * b) / 2 + \
              nm.abs(velo) * (a - b) / 2

        out[:] = 0.0

        # flux0 = (fl - fp)
        # flux1 = (- fl - fp + intg1)
        # out[:, 0, 0, 0] = -flux0
        # out[:, 0, 1, 0] = -flux1

        # for Legendre basis integral of higher order
        # functions of the basis is zero,
        # hence we calculate integral
        #
        # int_{j-1/2}^{j+1/2} f(u)dx
        #
        # only from the zero order function, over [-1, 1] - hence the 2
        intg1 = velo * u[:, 0] * 2
        intg2 = velo * u[:, 1] * 2 if n_el_nod > 2 else 0
        # i.e. intg1 = a * u0 * reference_el_vol

        flux = list()
        flux[:3] = (fl - fp), (- fl - fp + intg1), (fl - fp + intg2)
        for i in range(n_el_nod):
            out[:, :, i, 0] = -flux[i]


        status = None
        return status
import numpy as nm

from sfepy.discrete.dg.dg_basis import iter_by_order
from sfepy.discrete.dg.dg_field import get_raveler, get_unraveler
from sfepy.base.base import (get_default, output, assert_,
                             Struct, IndexedStruct)
from functools import  reduce
MACHINE_EPS = 1e-30


def minmod(a, b, c):
    """
    Minmod function of three variables, returns:
     _/ 0           , where sign(a) != sign(b) != sign(c)
      \ min(a,b,c)  , elsewhere
    :param a:
    :param b:
    :param c:
    :return:
    """
    seq = (nm.sign(a) == nm.sign(b)) & (nm.sign(b) == nm.sign(c))

    res = nm.zeros(nm.shape(a))
    res[seq] = nm.sign(a[seq]) * nm.minimum.reduce([nm.abs(b[seq]),
                                                    nm.abs(a[seq]),
                                                    nm.abs(c[seq])])

    return res

def minmod_seq(abc):
    seq =  nm.hstack([nm.sign(x) for x in abc])
    seq = seq[:, 0, None] == seq
    seq = seq.prod(axis=1).astype(bool)
    res = nm.zeros(nm.shape(abc[0]))
    res[seq] = nm.sign(abc[0][seq]) * nm.minimum.reduce([nm.abs(x[seq]) for x in abc])
    return res


class DGLimiter:
    name = "abstract DG limiter"

    def __init__(self, field, verbose=False):
        self.field = field
        self.n_el_nod = field.n_el_nod
        self.n_cell = field.n_cell
        self.ravel = get_raveler(self.n_el_nod, self.n_cell)
        self.unravel = get_unraveler(self.n_el_nod, self.n_cell)
        self.verbose = verbose

    def __call__(self, u):
        raise NotImplementedError("Called abstract limiter")


class IdentityLimiter(DGLimiter):
    name = "identity"

    def __call__(self, u):
        if self.verbose: output(self.name + " limiter")
        return u


class MomentLimiter1D(DGLimiter):
    """
    Krivodonova(2007): Limiters for high-order discontinuous Galerkin methods
    """
    name = "moment_1D_limiter"

    def __call__(self, u):
        """"
        :param u: solution at time step n in shape
        (order, n_space_nod)
        :return: limited solution
        """
        # for convenience do not try to limit FV
        if self.n_el_nod == 1:
            if self.verbose: output(self.name + " no limiting for FV.")
            return u
        u = self.unravel(u).swapaxes(0, 1)
        nbrhd_idx = self.field.get_facet_neighbor_idx()
        lidx = nbrhd_idx[:, 0, 0]
        ridx = nbrhd_idx[:, 1, 0]
        idx = nm.arange(nm.shape(u[0, 1:-1])[0])
        idx_bc = nm.arange(nm.shape(u[0, :])[0])

        nu = nm.copy(u)
        for l in range(self.n_el_nod - 1, 0, -1):
            tilu = minmod(nu[l, 1:-1][idx],
                          nu[l - 1, 2:][idx] - nu[l - 1, 1:-1][idx],
                          nu[l - 1, 1:-1][idx] - nu[l - 1, :-2][idx])

            # tilu_bc = minmod(nu[l, :][idx_bc],
            #               nu[l - 1, :][ridx] - nu[l - 1, :][idx_bc],
            #               nu[l - 1, :][idx_bc] - nu[l - 1, :][lidx])


            # idx_bc = nm.where(abs(tilu_bc - nu[l, :][idx_bc]) > MACHINE_EPS)[0]
            # lidx = lidx[idx_bc]
            # ridx = ridx[idx_bc]

            idx = nm.where(abs(tilu - nu[l, 1:-1][idx]) > MACHINE_EPS)[0]
            if self.verbose:
                output("{} limiting in {} cells out of {} :".
                       format(self.name, len(idx_bc), self.n_cell))
                output(idx_bc)
            if len(idx_bc) == 0:
                break
            nu[l, 1:-1][idx] = tilu[idx]
            # nu[l][idx_bc] = tilu_bc[idx_bc]

        return self.ravel(nu.swapaxes(0, 1))[:, 0]


class MomentLimiter2D(DGLimiter):
    """
    Krivodonova(2007): Limiters for high-order discontinuous Galerkin methods
    """
    name = "moment_limiter_2D"

    def __call__(self, u):
        if self.n_el_nod == 1:
            if self.verbose: output(self.name + " no limiting for FV.")
            return u


        nbrhd_idx = self.field.get_facet_neighbor_idx()
        inner_mask = nbrhd_idx[:, :,  0] > 0
        idx = nm.where(inner_mask.prod(axis=1))[0]
        nbrhd_idx = nbrhd_idx[:, :, 0]

        u = self.unravel(u).swapaxes(0, 1)
        nu = nm.zeros((self.field.approx_order + 1,) * 2 + u.shape[1:])
        tilu = nm.zeros(u.shape[1:])
        for l, (ii, jj) in enumerate(iter_by_order(self.field.approx_order, 2)):
            nu[ii, jj, ...] = u[l]

        for ii, jj in reversed(list(iter_by_order(self.field.approx_order, 2))):
            minmod_args = [nu[ii, jj, idx]]
            nbrhs = nbrhd_idx[idx]
            if ii - 1 >= 0:
                # right difference in x axis
                dx_r = nu[ii - 1, jj, nbrhs[:, 1]] - nu[ii - 1, jj, idx]
                # left differnce in x axis
                dx_l = nu[ii - 1, jj, idx] - nu[ii - 1, jj, nbrhs[:, 3]]
                minmod_args += [dx_r, dx_l]
            if jj - 1 >= 0:
                # right i.e. element "up" difference in y axis
                dy_up = nu[ii, jj - 1, nbrhs[:, 2]] - nu[ii, jj - 1 ,  idx]
                # left i.e. element "down" difference in y axis
                dy_dn = nu[ii, jj - 1,  idx] - nu[ii, jj - 1,  nbrhs[:, 0]]
                minmod_args += [dy_up, dy_dn]

            tilu[idx] = minmod_seq(minmod_args)
            idx = idx[nm.where(abs(tilu[idx] - nu[ii, jj, idx]) > MACHINE_EPS)[0]]
            if self.verbose:
                output("{} limiting {} in {} cells out of {} :".
                       format(self.name, (ii, jj) ,len(idx), self.n_cell))
                output(idx)
            if len(idx) == 0:
                break
            nu[ii, jj, idx] = tilu[idx]

        resu = nm.zeros(u.shape)
        for l, (ii, jj) in enumerate(iter_by_order(self.field.approx_order, 2)):
            resu[l] = nu[ii, jj]

        return self.ravel(resu.swapaxes(0, 1))[:, 0]



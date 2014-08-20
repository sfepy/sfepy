# -*- coding: utf-8 -*-
"""
Conversion of material parameters and other utilities.
"""
import os

import numpy as nm

from sfepy.base.base import Struct

def lame_from_youngpoisson(young, poisson, plane='strain'):
    r"""
    Compute Lamé parameters from Young's modulus and Poisson's ratio.

    The relationship between Lamé parameters and Young's modulus, Poisson's
    ratio (see [1],[2]):

    .. math::
        \lambda = {\nu E \over (1+\nu)(1-2\nu)},\qquad \mu = {E \over 2(1+\nu)}

    The plain stress hypothesis:

    .. math::
       \bar\lambda = {2\lambda\mu \over \lambda + 2\mu}

    [1] I.S. Sokolnikoff: Mathematical Theory of Elasticity. New York, 1956.

    [2] T.J.R. Hughes: The Finite Element Method, Linear Static and Dynamic
    Finite Element Analysis. New Jersey, 1987.
    """
    mu = young/(2.0*(1.0 + poisson))
    lam = young*poisson/((1.0 + poisson)*(1.0 - 2.0*poisson))

    if plane == 'stress':
        lam = 2*lam*mu/(lam + 2*mu)

    return lam, mu

def stiffness_from_lame(dim, lam, mu):
    r"""
    Compute stiffness tensor corresponding to Lamé parameters.

    .. math::
        {\bm D}_{(2D)} = \begin{bmatrix} \lambda + 2\mu & \lambda & 0\\
        \lambda & \lambda + 2\mu & 0\\ 0 & 0 & \mu \end{bmatrix}

    .. math::
        {\bm D}_{(3D)} = \begin{bmatrix} \lambda + 2\mu & \lambda &
        \lambda & 0 & 0 & 0\\ \lambda & \lambda + 2\mu & \lambda & 0 & 0 & 0 \\
        \lambda & \lambda & \lambda + 2\mu & 0 & 0 & 0 \\ 0 & 0 & 0 & \mu & 0 &
        0 \\ 0 & 0 & 0 & 0 & \mu & 0 \\ 0 & 0 & 0 & 0 & 0 & \mu\\ \end{bmatrix}
    """
    sym = (dim + 1) * dim / 2
    o = nm.array([1.] * dim + [0.] * (sym - dim), dtype=nm.float64)
    oot = nm.outer(o, o)[None, ...]
    do1 = nm.diag(o + 1.0)[None, ...]

    lam = nm.array(lam, ndmin=1)[:, None, None]
    mu = nm.array(mu, ndmin=1)[:, None, None]

    return (lam * oot + mu * do1).squeeze()

def stiffness_from_youngpoisson(dim, young, poisson, plane='strain'):
    """
    Compute stiffness tensor corresponding to Young's modulus and Poisson's
    ratio.
    """
    lam, mu = lame_from_youngpoisson(young, poisson, plane)

    return stiffness_from_lame(dim, lam, mu)

def stiffness_from_lame_mixed(dim, lam, mu):
    r"""
    Compute stiffness tensor corresponding to Lamé parameters for mixed
    formulation.

    .. math::
        {\bm D}_{(2D)} = \begin{bmatrix} \widetilde\lambda + 2\mu &
        \widetilde\lambda & 0\\ \widetilde\lambda & \widetilde\lambda + 2\mu &
        0\\ 0 & 0 & \mu \end{bmatrix}

    .. math::
        {\bm D}_{(3D)} = \begin{bmatrix} \widetilde\lambda + 2\mu &
        \widetilde\lambda & \widetilde\lambda & 0 & 0 & 0\\ \widetilde\lambda &
        \widetilde\lambda + 2\mu & \widetilde\lambda & 0 & 0 & 0 \\
        \widetilde\lambda & \widetilde\lambda & \widetilde\lambda + 2\mu & 0 &
        0 & 0 \\ 0 & 0 & 0 & \mu & 0 & 0 \\ 0 & 0 & 0 & 0 & \mu & 0 \\ 0 & 0 &
        0 & 0 & 0 & \mu\\ \end{bmatrix}

    where

    .. math::
       \widetilde\lambda = -{2\over 3} \mu
    """
    lam = - 2.0 / 3.0 * mu

    return stiffness_from_lame(dim, lam, mu)

def stiffness_from_youngpoisson_mixed(dim, young, poisson, plane='strain'):
    """
    Compute stiffness tensor corresponding to Young's modulus and Poisson's
    ratio for mixed formulation.
    """
    lam, mu = lame_from_youngpoisson(young, poisson, plane)

    return stiffness_from_lame_mixed(dim, lam, mu)

def bulk_from_lame(lam, mu):
    r"""
    Compute bulk modulus from Lamé parameters.

    .. math::
        \gamma = \lambda + {2 \over 3} \mu
    """
    return lam + 2.0 * mu / 3.0

def bulk_from_youngpoisson(young, poisson, plane='strain'):
    """
    Compute bulk modulus corresponding to Young's modulus and Poisson's ratio.
    """
    lam, mu = lame_from_youngpoisson(young, poisson, plane)

    return bulk_from_lame(lam, mu)

elastic_constants_relations = {
}

class ElasticConstants(Struct):
    r"""
    Conversion formulas for various groups of elastic constants. The elastic
    constants supported are:

      - :math:`E` : Young's modulus
      - :math:`\nu` : Poisson's ratio
      - :math:`K` : bulk modulus
      - :math:`\lambda` : Lamé's first parameter
      - :math:`\mu, G` : shear modulus, Lamé's second parameter
      - :math:`M` : P-wave modulus, longitudinal wave modulus

    The elastic constants are referred to by the following keyword arguments:
    young, poisson, bulk, lam, mu, p_wave.

    Exactly two of them must be provided to the __init__() method.

    Examples
    --------

     - basic usage::

        >>> from sfepy.mechanics.matcoefs import ElasticConstants
        >>> ec = ElasticConstants(lam=1.0, mu=1.5)
        >>> ec.young
        3.6000000000000001
        >>> ec.poisson
        0.20000000000000001
        >>> ec.bulk
        2.0
        >>> ec.p_wave
        4.0
        >>> ec.get(['bulk', 'lam', 'mu', 'young', 'poisson', 'p_wave'])
        [2.0, 1.0, 1.5, 3.6000000000000001, 0.20000000000000001, 4.0]

     - reinitialize existing instance::

        >>> ec.init(p_wave=4.0, bulk=2.0)
        >>> ec.get(['bulk', 'lam', 'mu', 'young', 'poisson', 'p_wave'])
        [2.0, 1.0, 1.5, 3.6000000000000001, 0.20000000000000001, 4.0]
    """
    def __init__(self, young=None, poisson=None, bulk=None, lam=None,
                 mu=None, p_wave=None, _regenerate_relations=False):
        """
        Set exactly two of the elastic constants, and compute the remaining.
        """
        self.names = ['bulk', 'lam', 'mu', 'young', 'poisson', 'p_wave']

        if _regenerate_relations:
            self.relations = self._construct_relations()

        else:
            import elastic_constants as ec
            self.relations = ec.relations
            self.ec = ec

        ## print sorted(self.relations.keys())
        ## print len(self.relations)

        self.init(young=young, poisson=poisson, bulk=bulk, lam=lam,
                  mu=mu, p_wave=p_wave)

    def _construct_relations(self):
        """
        Construct the dictionary of all relations among the six elastic
        constants and save it as `elastic_constants.py` module, that can be
        imported for reuse. Users should not call this!
        """
        import sympy as sm

        relations = {}

        def _expand_keys(sols):
            for key, val in sols.iteritems():
                if len(val) == 2 and (key.name == 'poisson'):
                    val = val[0]
                else:
                    val = val[-1]
                skey = tuple(sorted([ii.name for ii in val.atoms()
                                     if ii.is_Symbol])) + (key.name,)
                if skey in relations:
                    print '!', skey
                relations[skey] = val

        bulk, lam, mu, young, poisson, p_wave = sm.symbols(self.names,
                                                           real=True)

        _expand_keys(sm.solve(bulk - (lam + 2 * mu / 3)))
        _expand_keys(sm.solve(young - (mu * (3 * lam + 2 * mu) / (lam + mu))))
        _expand_keys(sm.solve(poisson - (lam / (2 * (lam + mu)))))
        _expand_keys(sm.solve(p_wave - (lam + 2 * mu)))

        _expand_keys(sm.solve(bulk - (young / (3 * (1 - 2 * poisson)))))
        _expand_keys(sm.solve(p_wave - ((young * (1 - poisson))
                                        / ((1 + poisson) * (1 - 2 * poisson)))))
        # Choose the correct root manually.
        ## relations[('p_wave', 'young', 'poisson')] \
        ##                 = (young - p_wave + (-10*p_wave*young + young**2 +
        ##                                      9*p_wave**2)**(0.5))/(4*p_wave)
        _expand_keys(sm.solve(lam - (young * poisson
                                     / ((1 + poisson) * (1 - 2 * poisson)))))
        # Choose the correct root.
        ## relations[('lam', 'young', 'poisson')] \
        ##                   = (lam + young - (2*lam*young + young**2 +
        ##                                     9*(lam**2))**(0.5))/(-4*lam)
        _expand_keys(sm.solve(mu - (young / (2 * (1 + poisson)))))

        _expand_keys(sm.solve(bulk - (young * mu / (3 * (3 * mu - young)))))
        _expand_keys(sm.solve(p_wave - (mu * (4 * mu - young)
                                        / (3 * mu - young))))

        _expand_keys(sm.solve(young - (9 * bulk * (bulk - lam)
                                       / (3 * bulk - lam))))
        _expand_keys(sm.solve(poisson - (lam / (3 * bulk - lam))))
        _expand_keys(sm.solve(p_wave - (3 * bulk - 2 * lam)))

        _expand_keys(sm.solve(poisson - ((3 * bulk - 2 * mu)
                                         / (2 * (3 * bulk + mu)))))
        _expand_keys(sm.solve(p_wave - (bulk + 4 * mu / 3)))

        _expand_keys(sm.solve(p_wave - (lam * (1 - poisson) / poisson)))

        _expand_keys(sm.solve(p_wave - (2 * mu * (1 - poisson)
                                        / (1 - 2 * poisson))))

        _expand_keys(sm.solve(p_wave - (3 * bulk * (1 - poisson)
                                        / (1 + poisson))))

        _expand_keys(sm.solve(p_wave - (3 * bulk * (3 * bulk + young)
                                        / (9 * bulk - young))))

        _expand_keys(sm.solve(young - ((lam*p_wave + p_wave**2 - 2*lam**2)
                                       / (lam + p_wave))))

        fd = open(os.path.join(os.path.dirname(__file__),
                               'elastic_constants.py'), 'w')
        fd.write("""
from __future__ import division

import sympy as sm

names = ['bulk', 'lam', 'mu', 'young', 'poisson', 'p_wave']
bulk, lam, mu, young, poisson, p_wave = sm.symbols(names, real=True)

relations = {
%s
}
        """ % ',\n'.join(['    %s : %s' % (key, val)
                         for key, val in relations.iteritems()]))
        fd.close()

        return relations

    def init(self, young=None, poisson=None, bulk=None, lam=None,
             mu=None, p_wave=None):
        """
        Set exactly two of the elastic constants, and compute the
        remaining. (Re)-initializes the existing instance of ElasticConstants.
        """
        Struct.__init__(self, young=young, poisson=poisson, bulk=bulk, lam=lam,
                        mu=mu, p_wave=p_wave)

        values = {}
        for key, val in self.__dict__.iteritems():
            if (key in self.names) and (val is not None):
                sym = getattr(self.ec, key)
                values[sym] = val

        known = values.keys()
        if len(known) != 2:
            raise ValueError('exactly two elastic constants must be provided!')
        known = [ii.name for ii in known]

        unknown = set(self.names).difference(known)

        for name in unknown:
            key = tuple(sorted(known)) + (name,)
            val = float(self.relations[key].n(subs=values))
            setattr(self, name, val)

    def get(self, names):
        """
        Get the named elastic constants.
        """
        out = [getattr(self, name) for name in names]
        return out

class TransformToPlane(Struct):
    """
    Transformations of constitutive law coefficients of 3D problems to 2D.
    """

    def __init__(self, iplane=None):
        """
        Parameters
        ----------
        iplane : list
            The vector of indices denoting the plane, e.g.: [0, 1]
        """
        if iplane is None:
            iplane = [0, 1]

        # Choose the "master" variables and the "slave" ones
        # ... for vectors
        i_m = nm.sort(iplane)
        i_s = nm.setdiff1d(nm.arange(3), i_m)

        # ... for second order tensors (symmetric storage)
        i_ms = {(0, 1) : [0, 1, 3],
                (0, 2) : [0, 2, 4],
                (1, 2) : [1, 2, 5]}[tuple(i_m)]
        i_ss = nm.setdiff1d(nm.arange(6), i_ms)

        Struct.__init__(self, iplane=iplane,
                        i_m=i_m, i_s=i_s, i_ms=i_ms, i_ss=i_ss)

    def tensor_plane_stress(self, c3=None, d3=None, b3=None):
        """
        Transforms all coefficients of the piezoelectric constitutive law
        from 3D to plane stress problem in 2D: strain/stress ordering: 11 22
        33 12 13 23. If `d3` is None, uses only the stiffness tensor `c3`.

        Parameters
        ----------
        c3 : array
            The stiffness tensor.
        d3 : array
            The dielectric tensor.
        b3 : array
            The piezoelectric coupling tensor.
        """
        mg = nm.meshgrid

        cs = c3[mg(self.i_ss, self.i_ss)]
        cm = c3[mg(self.i_ss, self.i_ms)].T
        if d3 is None: # elasticity only.
            A = cs
            Feps = cm

            Ainv = nm.linalg.inv(A)
            c2 = c3[mg(self.i_ms, self.i_ms)] \
                 - nm.dot(Feps.T, nm.dot(Ainv, Feps))

            return c2

        else:
            dm = d3[mg(self.i_s, self.i_m)].T
            ds = d3[mg(self.i_s, self.i_s)]

            ii = mg(self.i_s, self.i_ss)
            A = nm.r_[nm.c_[cs, b3[ii]],
                      nm.c_[b3[ii].T, -ds]] #=> sym !!!
            F = nm.r_[nm.c_[cm, b3[mg(self.i_m, self.i_ss)]],
                      nm.c_[b3[mg(self.i_s, self.i_ms)].T, -dm]]

            Feps = F[:, :3]
            FE = F[:, 3:]
            Ainv = nm.linalg.inv(A)
            c2 = c3[mg(self.i_ms, self.i_ms)] \
                 - nm.dot(Feps.T, nm.dot(Ainv, Feps))
            d2 = d3[mg(self.i_m, self.i_m)] \
                 - nm.dot(FE.T, nm.dot(Ainv, FE))
            b2 = b3[mg(self.i_m, self.i_ms)].T \
                 - nm.dot(FE.T, nm.dot(Ainv, Feps))

            return c2, d2, b2

class ElasticTensor():
    """
    This class represents linear fourth-order tensor of elastic parameters
    for both stiffness and compliance. It also evaluates the tensor as a matrix
    using engineering Mandel's or Voight's notation. The plane stress or strain
    is also available.
    """

    def __init__(self, bulk=None, mu=None, stiffness=True, plane=None):
        """
        Parameters
        ----------
        bulk : float
            bulk modulus
        mu : float
            shear modulus
        stiffness : boolean
            determines whether the values are evaluated for stiffness (True)
            or for its inverse compliance (False)
        plane : None, 'strain', or 'stress'
            determines the physical problem; None represents three-dimensional
            problem and 'strain' or 'stress' represent reduced problem of plane
            strain or stress respectively
        """
        if stiffness:
            self.stiffness = stiffness
            self.val_type = 'stiffness'
        else:
            self.val_type = 'compliance'

        if plane is None:
            self.dim = 3
        elif plane in ['stress', 'strain']:
            self.dim = 2
        else:
            raise ValueError("This type of plane (%s) is not supported." \
                            % str(plane))

        self.plane = plane

        self.sym = self.get_sym(self.dim)
        self.bulk = bulk
        self.mu = mu

        # set up values to tensor and matrices
        _, volumetric, deviatoric = self.get_decomposition()
        if stiffness:
            self.val = 3.*bulk*volumetric + 2.*mu*deviatoric
        else:
            self.val = 1./(3*bulk)*volumetric + 1./(2*mu)*deviatoric

        self.mandel = self.create_mandel(self.val)
        self.voight = self.create_voight(self.val)

        if plane is not None:
            self.val = self.val[0:2, 0:2, 0:2, 0:2]

        if (stiffness and plane == 'strain') or \
                (not stiffness and plane == 'stress'):
            self.mandel = self.get_plane(self.mandel)
            self.voight = self.get_plane(self.voight)
        elif (not stiffness and plane == 'strain') or \
                (stiffness and plane == 'stress'):
            inv = nm.linalg.inv
            self.mandel = inv(self.get_plane(inv(self.mandel)))
            self.voight = inv(self.get_plane(inv(self.voight)))
        else:
            pass

    @staticmethod
    def get_plane(val, ind=None):
        if ind is None:
            ind = [0, 1]
        ind_shear = range(3)
        ind_shear.remove(ind[0])
        ind_shear.remove(ind[1])
        ind.append(ind_shear[0]+3)
        mat = val[ind][:, ind]
        return mat

    def __repr__(self):
        ss = "Class: %s\n" % (self.__class__.__name__)
        ss += '    stiffness = %s (%s)\n' % (self.stiffness, self.val_type)
        ss += '    dim = %d' % (self.dim)
        if self.plane is None:
            ss += '\n'
        else:
            ss += ' (plane %s)\n' % self.plane
        ss += '    bulk = %s\n' % str(self.bulk)
        ss += '    mu = %s\n' % str(self.mu)
        return ss

    @staticmethod
    def get_sym(dim):
        return dim*(dim+1)/2

    @staticmethod
    def get_decomposition():
        """
        It produces symmetrized fourth-order identity, hydrostatic, and
        deviatoric projections.

        Returns
        -------
        idsym : numpy.array of shape = (3, 3, 3, 3)
            symmetrized identity operator
        volumetric : numpy.array of shape = (3, 3, 3, 3)
            hydrostatic projection
        deviatoric : numpy.array of shape = (3, 3, 3, 3)
            deviatoric projection
        """
        id = nm.eye(3)
        volumetric = 1./3*nm.einsum('ij,kl', id, id)
        idsym = 0.5*(nm.einsum('ik,jl', id, id) + nm.einsum('il,jk', id, id))
        deviatoric = idsym - volumetric
        return idsym, volumetric, deviatoric

    @staticmethod
    def create_mandel(mat):
        """
        It transfer symmetric four-order tensor (or matrix) to
        second order tensor (or vector) using Mandel's notation.

        Parameters
        ----------
        mat : numpy.array of shape = (d, d, d, d) or (d, d) for dimension d
            fourth-order tensor of elastic parameters

        Returns
        -------
        vec : numpy.array of shape = (sym, sym) or (sym,) for sym = d*(d+1)/2
            second-order tensor of elastic parameters with Mandel's notation
        """
        dim = mat.shape[0]
        sym = dim*(dim+1)/2
        if mat.ndim == 4:
            vec = nm.zeros([sym, sym], dtype=mat.dtype)
            for ii in nm.arange(dim):
                for jj in nm.arange(dim):
                    kk = range(dim)
                    kk.remove(ii)
                    ll = range(dim)
                    ll.remove(jj)
                    vec[ii, jj] = mat[ii, ii, jj, jj]
                    vec[ii, jj+dim] = 2**.5*mat[ii, ii, ll[0], ll[1]]
                    vec[jj+dim, ii] = vec[ii, jj+dim]
                    vec[ii+dim, jj+dim] = 2*mat[kk[0], kk[1], ll[0], ll[1]]
        elif mat.ndim == 2:
            vec = nm.zeros(sym, dtype=mat.dtype)
            vec[:dim] = nm.diag(mat)
            if dim == 2:
                vec[dim] = 2**.5*mat[0, 1]
            elif dim == 3:
                for ii in nm.arange(sym-dim):
                    ind = range(sym-dim)
                    ind.remove(ii)
                    vec[dim+ii] = 2**.5*mat[ind[0], ind[1]]
            else:
                raise ValueError("Incorrect dimension (%d)" % dim)
        return vec

    @staticmethod
    def dispose_mandel(vec):
        vec = vec.squeeze()
        sym = vec.shape[0]
        dimfun = lambda sym: int((-1.+(1+8*sym)**.5)/2)
        dim = dimfun(sym)

        if vec.ndim == 2: # matrix
            assert_(vec.shape[0]==vec.shape[1])
            assert_(dim==3, 'test')
            mat = nm.zeros(dim*nm.ones(2*vec.ndim))

            for ii in nm.arange(dim):
                for jj in nm.arange(dim):
                    kk = range(dim)
                    kk.remove(ii)
                    ll = range(dim)
                    ll.remove(jj)
                    mat[ii, ii, jj, jj] = vec[ii, jj]
                    mat[ii, ii, ll[0], ll[1]] = vec[ii, jj+dim] / 2**.5
                    mat[ll[0], ll[1], ii, ii] = mat[ii, ii, ll[0], ll[1]]
                    mat[kk[0], kk[1], ll[0], ll[1]] = vec[ii+dim, jj+dim] / 2.
                    mat[kk[1], kk[0], ll[0], ll[1]] = vec[ii+dim, jj+dim] / 2.
                    mat[kk[0], kk[1], ll[1], ll[0]] = vec[ii+dim, jj+dim] / 2.
                    mat[kk[1], kk[0], ll[1], ll[0]] = vec[ii+dim, jj+dim] / 2.

        elif vec.ndim == 1: # vector
            mat = nm.diag(vec[:dim])

            if dim == 2:
                mat[0, 1] = vec[-1]/2**0.5
                mat[1, 0] = vec[-1]/2**0.5
            elif dim == 3:
                for ii in nm.arange(sym-dim):
                    ind = range(sym-dim)
                    ind.remove(ii)
                    mat[ind[0], ind[1]] = vec[dim+ii]/2.**.5
                    mat[ind[1], ind[0]] = vec[dim+ii]/2.**.5
            else:
                raise ValueError("Incorrect dimension (%d)" % dim)

        return mat

    @staticmethod
    def create_voight(mat, valtype='strain'):
        """
        It transfer symmetric four-order tensor to second order tensor
        using Voight's notation.

        Parameters
        ----------
        mat : numpy.array of shape = (3, 3, 3, 3)
            fourth-order tensor of elastic parameters
        valtype : one of 'strain' or 'stress'
            this distinguish a engineering notation for strain and stress

        Returns
        -------
        vec : numpy.array of shape = (6, 6)
            second-order tensor of elastic parameters with Voight's notation
        """
        dim = mat.shape[0]
        sym = dim*(dim+1)/2
        if mat.ndim == 4:
            vec = nm.zeros([sym, sym], dtype=mat.dtype)
            for ii in nm.arange(dim):
                for jj in nm.arange(dim):
                    kk = range(dim)
                    kk.remove(ii)
                    ll = range(dim)
                    ll.remove(jj)
                    vec[ii, jj] = mat[ii, ii, jj, jj]
                    vec[ii, jj+dim] = mat[ii, ii, ll[0], ll[1]]
                    vec[jj+dim, ii] = vec[ii, jj+dim]
                    vec[ii+dim, jj+dim] = mat[kk[0], kk[1], ll[0], ll[1]]

        elif mat.ndim == 2:
            vec = nm.zeros(sym, dtype=mat.dtype)
            vec[:dim] = nm.diag(mat)
            if valtype == 'strain':
                coef = 2.
            elif valtype == 'stress':
                coef = 1.
            else:
                msg = "Parameter valtype (%s) should be one of 'strain' or\
                     'stress'." % (str(valtype),)
                raise ValueError(msg)

            if dim == 2:
                vec[dim] = coef*mat[0, 1]
            elif dim == 3:
                for ii in nm.arange(sym-dim):
                    ind = range(sym-dim)
                    ind.remove(ii)
                    vec[dim+ii] = coef*mat[ind[0], ind[1]]
            else:
                raise ValueError("Incorrect dimension (%d)" % dim)
        return vec
>>>>>>> new general term

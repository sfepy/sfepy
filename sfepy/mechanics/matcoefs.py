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
    sym = (dim + 1) * dim // 2
    o = nm.array([1.] * dim + [0.] * (sym - dim), dtype=nm.float64)
    oot = nm.outer(o, o)
    do1 = nm.diag(o + 1.0)

    lam = nm.array(lam)[..., None, None]
    mu = nm.array(mu)[..., None, None]
    return (lam * oot + mu * do1)

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

def lame_from_stiffness(stiffness, plane='strain'):
    """
    Compute Lamé parameters from an isotropic stiffness tensor.
    """
    lam = stiffness[..., 0, 1]
    mu = stiffness[..., -1, -1]
    if plane == 'stress':
        lam = - 2.0 * mu * lam / (lam - 2.0 * mu)

    return lam, mu

def youngpoisson_from_stiffness(stiffness, plane='strain'):
    """
    Compute Young's modulus and Poisson's ratio from an isotropic stiffness
    tensor.
    """
    lam, mu = lame_from_stiffness(stiffness, plane=plane)
    young = (3*lam*mu + 2*mu**2) / (lam + mu)
    poisson = lam / (2*lam + 2*mu)

    return young, poisson


def stiffness_from_yps_ortho3(young, poisson, shear):
    r"""
    Compute 3D stiffness tensor :math:`{\bm D}` of an orthotropic linear
    elastic material. Young's modulus (:math:`[E_1, E_2, E_3]`),
    Poisson's ratio (:math:`[\nu_{12}, \nu_{13}, \nu_{23}]`),
    and shear modulus (:math:`[G_{12}, G_{13}, G_{23}]`) are given.

    .. math::
        {\bm C}_{(3D)} = \begin{bmatrix}
        1/E_1 & -\nu_{21}/E_2 & -\nu_{31}/E_3 & 0 & 0 & 0\\
        -\nu_{12}/E_1 & 1/E_2 & -\nu_{32}/E_3 & 0 & 0 & 0\\
        -\nu_{13}/E_1 & -\nu_{23}/E_2 & 1/E_3 & 0 & 0 & 0\\
        0 & 0 & 0 & 1/G_{12} & 0 & 0 \\
        0 & 0 & 0 & 0 & 1/G_{13} & 0 \\
        0 & 0 & 0 & 0 & 0 & 1/G_{23} \end{bmatrix}

    .. math::
        {\bm D}_{(3D)} = \mathrm{inv}({\bm C}_{(3D)})

    .. math::
        \nu_{21} = \nu_{12}\frac{E_2}{E_1},\quad
        \nu_{31} = \nu_{13}\frac{E_3}{E_1},\quad
        \nu_{32} = \nu_{23}\frac{E_3}{E_2}

    [1] R.M. Jones: Mechanics of composite materials. 1999.
    """
    # symmetric order [11, 22, 33, 12, 13, 23]
    compl = nm.zeros((6, 6), dtype=nm.float64)

    compl[0, 0] = 1. / young[0]
    compl[1, 1] = 1. / young[1]
    compl[2, 2] = 1. / young[2]
    compl[0, 1] = compl[1, 0] = -poisson[0] / young[0]
    compl[0, 2] = compl[2, 0] = -poisson[1] / young[0]
    compl[1, 2] = compl[2, 1] = -poisson[2] / young[1]
    compl[3, 3] = 1. / shear[0]
    compl[4, 4] = 1. / shear[1]
    compl[5, 5] = 1. / shear[2]

    stiff = nm.linalg.inv(compl)

    return stiff


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
            from . import elastic_constants as ec
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
            for key, val in sols.items():
                if len(val) == 2 and (key.name == 'poisson'):
                    val = val[0]
                else:
                    val = val[-1]
                skey = tuple(sorted([ii.name for ii in val.atoms()
                                     if ii.is_Symbol])) + (key.name,)
                if skey in relations:
                    print('!', skey)
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
import sympy as sm

names = ['bulk', 'lam', 'mu', 'young', 'poisson', 'p_wave']
bulk, lam, mu, young, poisson, p_wave = sm.symbols(names, real=True)

relations = {
%s
}
        """ % ',\n'.join(['    %s : %s' % (key, val)
                         for key, val in relations.items()]))
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
        for key, val in self.__dict__.items():
            if (key in self.names) and (val is not None):
                sym = getattr(self.ec, key)
                values[sym] = val

        known = list(values.keys())
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

        cs = c3[tuple(mg(self.i_ss, self.i_ss))]
        cm = c3[tuple(mg(self.i_ss, self.i_ms))].T
        if d3 is None: # elasticity only.
            A = cs
            Feps = cm

            Ainv = nm.linalg.inv(A)
            c2 = c3[tuple(mg(self.i_ms, self.i_ms))] \
                 - nm.dot(Feps.T, nm.dot(Ainv, Feps))

            return c2

        else:
            dm = d3[tuple(mg(self.i_s, self.i_m))].T
            ds = d3[tuple(mg(self.i_s, self.i_s))]

            ii = tuple(mg(self.i_s, self.i_ss))
            A = nm.r_[nm.c_[cs, b3[ii]],
                      nm.c_[b3[ii].T, -ds]] #=> sym !!!
            F = nm.r_[nm.c_[cm, b3[tuple(mg(self.i_m, self.i_ss))]],
                      nm.c_[b3[tuple(mg(self.i_s, self.i_ms))].T, -dm]]

            Feps = F[:, :3]
            FE = F[:, 3:]
            Ainv = nm.linalg.inv(A)
            c2 = c3[tuple(mg(self.i_ms, self.i_ms))] \
                 - nm.dot(Feps.T, nm.dot(Ainv, Feps))
            d2 = d3[tuple(mg(self.i_m, self.i_m))] \
                 - nm.dot(FE.T, nm.dot(Ainv, FE))
            b2 = b3[tuple(mg(self.i_m, self.i_ms))].T \
                 - nm.dot(FE.T, nm.dot(Ainv, Feps))

            return c2, d2, b2

def youngpoisson_from_wave_speeds(vp, vs, rho):
    r"""
    Compute the Young's modulus :math:`E` and Poisson's ratio :math:`\nu` from
    the P- and S-wave speeds in a homogeneous isotropic material.

    .. math::
       E = {\rho v_s^2 (3 v_p^2 - 4 v_s^2) \over (v_p^2 - v_s^2)}

    .. math::
       \nu = {(v_p^2/2 - v_s^2) \over (v_p^2 - v_s^2)}

    Parameters
    ----------
    vp : float or array
        The P-wave speed.
    vs : float or array
        The S-wave speed.
    rho : float or array
        The density.

    Returns
    -------
    young : float or array
        The Young's modulus.
    poisson : float or array
        The Poisson's ratio.
    """
    vs2 = vs**2
    vp2 = vp**2
    aux = vp2 - vs2
    return (vs2 * rho * (3.0*vp2 - 4.0*vs**2)/aux,
            (vp2/2.0 - vs2)/aux)

def wave_speeds_from_youngpoisson(young, poisson, rho):
    r"""
    Compute the P- and S-wave speeds from the Young's modulus :math:`E` and
    Poisson's ratio :math:`\nu` in a homogeneous isotropic material.

    .. math::
       v_p^2 = {E (1 - \nu) \over \rho (1 + \nu) (1 - 2 \nu)}
             = {(\lambda + 2 \mu) \over \rho}

    .. math::
       v_s^2 = {E \over 2 \rho (1 + \nu)}
             = {\mu \over \rho}

    Parameters
    ----------
    young : float or array
        The Young's modulus.
    poisson : float or array
        The Poisson's ratio.
    rho : float or array
        The density.

    Returns
    -------
    vp : float or array
        The P-wave speed.
    vs : float or array
        The S-wave speed.
    """
    return (nm.sqrt(young * (1.0 - poisson)
                    / (rho * (1.0 + poisson) * (1.0 - 2.0*poisson))),
            nm.sqrt(young / (2.0 * rho * (1.0 + poisson))))

# -*- coding: utf-8 -*-
"""
Terms implementing shell elements.
"""
import numpy as nm

from sfepy.terms.terms import Term
from sfepy.linalg import dot_sequences as ddot
import sfepy.mechanics.shell10x as shell10x

class Shell10XTerm(Term):
    r"""
    The shell10x element term based on the Reissner-Mindlin theory [1]_, [2]_,
    corresponding to a shell of thickness :math:`t`.

    The term requires a custom 3D quadrature, where the :math:`z` components of
    quadrature point coordinates are transformed from :math:`[0, 1]` to
    :math:`[-t/2, t/2]`, and the quadrature weights are multiplied by
    :math:`t`. The variables :math:`\ul{v}` and :math:`\ul{u}` have to use
    :class:`Shell10XField <sfepy.discrete.structural.fields.Shell10XField>` and
    have six components. The reference element mapping is implemented by
    :class:`Shell10XMapping
    <sfepy.discrete.structural.mappings.Shell10XMapping>`. The term does not
    implement the piezo-electric components of the shell10x element yet.

    The term has to be used with quadrilateral cells in 3D and should behave as
    the linear elastic term, but with fewer degrees of freedom for the same
    accuracy for shell-like structures. The shell has six degrees of freedom in
    each of the four nodes: :math:`{\mathbf u}_i = [ u_i, v_i, w_i, \alpha_i,
    \beta_i, \gamma_i ]^T`, :math:`i = 1, 2, 3, 4`. The strain and stress
    vectors are calculated in a local (co-rotational) coordinate system given
    by basis vectors :math:`{\mathbf e}'_1`, :math:`{\mathbf e}'_2` and
    :math:`{\mathbf e}'_3`. It holds that

    .. math::

       [ u'_i, v'_i, w'_i, \alpha'_i, \beta'_i, \gamma'_i ]^T
       = {\hat{\mathbf H}}^T{\mathbf u}_i\,

    where

    .. math::

       {\hat{\mathbf H}} = \left[ \begin{array}{cc} {\mathbf H} & \\ & {\mathbf
       H} \\ \end{array} \right]\qquad
       {\mathrm{and}} \qquad {\mathbf H} =
       [{\mathbf e}'_1 \, {\mathbf e}'_2 \, {\mathbf e}'_3]

    is a nodal DOF transformation matrix.

    The local displacements :math:`u'`, :math:`v'` and :math:`w'` at any point
    in the layer characterized by the isoparametric coordinates :math:`\xi`,
    :math:`\eta` and :math:`\zeta` (:math:`\xi, \eta, \zeta \in
    \left\langle-1,1\right\rangle`) are interpolated from the nodal
    displacement and rotation values (i.e. both membrane and bending
    components) using standard isoparametric approximation functions for a
    quadrilateral, hence

    .. math::

       \begin{array}{lcl}
       u'(\xi,\eta,\zeta) &=& \sum\limits_{i=1}^{4}
       N_i(\xi,\eta)\cdot ( u'_i + \bar{u}_i )\, , \\
       v'(\xi,\eta,\zeta) &=&
       \sum\limits_{i=1}^{4} N_i(\xi,\eta)\cdot ( v'_i + \bar{v}_i )\, , \\
       w'(\xi,\eta,\zeta) &=& \sum\limits_{i=1}^{4} N_i(\xi,\eta)\cdot
       ( w'_i + \bar{w}_i )\,
       \end{array}

    where :math:`\tilde{u}_i`, :math:`\tilde{v}_i` and :math:`\tilde{w}_i` are
    the bending components of displacements calculated from displacements due
    to rotations :math:`\tilde{\alpha}_i` and :math:`\tilde{\beta}_i` about
    local nodal axes :math:`\tilde{\mathbf{e}}_i` as

    .. math::

       \left[ \begin{array}{c} \bar{u}\\ \bar{v}\\ \bar{w}\\ \end{array}
       \right]_i = \tilde{\zeta} \left[ \begin{array}{cc} & \\
       \tilde{\mathbf{e}}_1 & -\tilde{\mathbf{e}}_2 \\ & \\ \end{array}
       \right]_i \left[ \begin{array}{c}
       \tilde{\mathbf{e}}_2^T \\ \tilde{\mathbf{e}}_1^T \\ \end{array}
       \right]_i \left[ \begin{array}{c} \alpha' \\ \beta' \\ \gamma' \\
       \end{array} \right]_i \,

    where :math:`\tilde{\zeta} = (t/2)\zeta`. The local nodal axes
    :math:`\tilde{\mathbf{e}}_i` are constructed in order to describe the
    behavior of warped (non-planar) elements adequately.

    The term employs three shell element enhancements:

    - DSG method
    - EAS method
    - drilling rotations lock (parameter :math:`\chi` - a good value is about
      :math:`10^{-7}`)

    For detailed theoretical information see the references.

    .. [1] Zemčík, R., Rolfes, R., Rose, M. and Tessmer, J. (2006),
    High-Performance 4-Node Shell Element with Piezoelectric Coupling Mechanics
    of Advanced Materials and Structures Vol. 13, Iss. 5,
    doi:10.1080/15376490600777657

    .. [2] Zemčík, R., Rolfes, R., Rose, M. and Teßmer, J. (2007),
    High-performance four-node shell element with piezoelectric coupling for
    the analysis of smart laminated structures. Int. J. Numer. Meth. Engng.,
    70: 934–961. doi:10.1002/nme.1909

    :Definition:

    .. math::
        \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})

    :Arguments:
        - material_d     : :math:`D`
        - material_drill : :math:`\chi`
        - virtual        : :math:`\ul{v}`
        - state          : :math:`\ul{u}`
    """
    name = 'dw_shell10x'
    arg_types = ('material_d', 'material_drill', 'virtual', 'state')
    arg_shapes = {'material_d' : '6, 6', 'material_drill' : '.: 1',
                  'virtual' : (6, 'state'), 'state' : 6}
    geometries = ['3_2_4']
    integration = 'custom'
    poly_space_basis = 'shell10x'

    def set_integral(self, integral):
        """
        Set the term integral.
        """
        if (integral.mode != 'custom') or (integral.coors.shape[1] != 3):
            raise ValueError("dw_shell10x term requires 'custom' 3D integral!")

        Term.set_integral(self, integral)

    def get_physical_qps(self):
        """
        Get physical quadrature points corresponding to the term region
        and integral.
        """
        mapping, _ = self.get_mapping(self.get_virtual_variable())

        qp_coors, _ = self.integral.get_qp('')
        phys_qps = mapping.get_physical_qps(qp_coors)

        return phys_qps

    @staticmethod
    def function(out, mtx_k, el_u, fmode):
        if fmode == 0:
            out[:, 0, ...] = ddot(mtx_k, el_u[..., None], 'AB')

        else:
            out[:, 0, ...] = mtx_k

        return 0

    def get_fargs(self, mtx_d, drill, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        from sfepy.discrete.variables import create_adof_conn

        vv, vu = virtual, state
        geo, _ = self.get_mapping(vv)

        # Displacements of element nodes.
        vec_u = vu()
        econn = vu.field.get_econn('cell', self.region)
        adc = create_adof_conn(nm.arange(vu.n_dof, dtype=nm.int32), econn, 6, 0)
        el_u = vec_u[adc]

        qp_coors = geo.qp.vals

        dsg = shell10x.get_dsg_strain(geo.coors_loc, geo.qp.vals)
        mtx_b = shell10x.create_strain_matrix(geo.bfgm, geo.dxidx, dsg)
        mtx_be = shell10x.add_eas_dofs(mtx_b, qp_coors, geo.det,
                                       geo.det0, geo.dxidx0)

        mtx_dr = shell10x.rotate_elastic_tensor(mtx_d, geo.bfu, geo.ebs)
        mtx_k_qp = ddot(ddot(mtx_be, mtx_dr, 'ATB'), mtx_be, 'AB')

        # Integrate.
        mtx_k = (mtx_k_qp * (geo.det * geo.qp.weights)[..., None, None]).sum(1)

        # Static condensation.
        k11 = mtx_k[:, :24, :24]
        k12 = mtx_k[:, :24, 24:]
        k22 = mtx_k[:, 24:, 24:]

        mtx_k = k11 - ddot(ddot(k12, nm.linalg.inv(k22)),
                           k12.transpose(0, 2, 1))

        if drill != 0.0:
            coefs =  mtx_dr[..., 3, 3].mean(1) * geo.volume * drill
            mtx_k = shell10x.lock_drilling_rotations(mtx_k, geo.ebs, coefs)

        # DOFs in mtx_k are DOF-by-DOF. Transform is u and phi node-by-node.
        blocks = nm.arange(24).reshape(2, 3, 4)
        blocks = blocks.transpose((0, 2, 1)).reshape((-1, 3))
        shell10x.transform_asm_matrices(mtx_k[:, None, ...], geo.mtx_t, blocks)

        fmode = diff_var is not None

        return (mtx_k, el_u, fmode)

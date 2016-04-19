"""
Terms implementing shell elements.
"""
import numpy as nm

from sfepy.terms.terms import Term
from sfepy.linalg import dot_sequences as ddot
import sfepy.mechanics.shell10x as shell10x

class Shell10XTerm(Term):
    r"""
    The shell10x element term.

    :Arguments:
        - material_d     : :math:`D`
        - material_drill : :math:`\chi`
        - virtual        : :math:`v`
        - state          : :math:`u`
    """
    name = 'dw_shell10x'
    arg_types = ('material_d', 'material_drill', 'virtual', 'state')
    arg_shapes = {'material_d' : '6, 6', 'material_drill' : '.: 1',
                  'virtual' : (6, 'state'), 'state' : 6}
    geometries = ['3_2_4']
    integration = 'custom'
    poly_space_base = 'shell10x'

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
        econn = vu.field.get_econn('volume', self.region)
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

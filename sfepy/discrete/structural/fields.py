"""
Fields corresponding to structural elements.
"""
from sfepy.base.base import Struct
from sfepy.discrete.fem.fields_nodal import H1NodalMixin, VolumeField
from sfepy.discrete.structural.mappings import Shell10XMapping
from sfepy.discrete import PolySpace

class Shell10XField(H1NodalMixin, VolumeField):
    """
    The field for the shell10x element.
    """
    family_name = 'volume_H1_shell10x'

    def _create_interpolant(self):
        name = '%s_%s_%s_%d' % (self.gel.name, self.space,
                                  self.poly_space_base, self.approx_order)
        ps = PolySpace.any_from_args(name, self.gel, self.approx_order,
                                     base='lagrange', force_bubble=False)
        self.poly_space = ps

    def create_mapping(self, region, integral, integration,
                       return_mapping=True):
        """
        Create a new reference mapping.
        """
        if integration != 'custom':
            msg = "Shell10XField requires 'custom' integration instead of '%s'!"
            raise ValueError(msg % integration)

        qp = self.get_qp('v', integral)

        mapping = Shell10XMapping(region, self)
        out = mapping.get_mapping(qp.vals, qp.weights)

        # Store the integral used.
        out.integral = integral
        out.qp = qp

        if return_mapping:
            out = (out, mapping)

        return out

    def create_output(self, dofs, var_name, dof_names=None,
                      key=None, thickness=None, **kwargs):
        """
        Convert the DOFs corresponding to the field to a dictionary of
        output data usable by Mesh.write().

        Parameters
        ----------
        dofs : array, shape (n_nod, n_component)
            The array of DOFs reshaped so that each column corresponds
            to one component.
        var_name : str
            The variable name corresponding to `dofs`.
        dof_names : tuple of str
            The names of DOF components.
        key : str, optional
            The key to be used in the output dictionary instead of the
            variable name.

        Returns
        -------
        out : dict
            The output dictionary.
        """
        out = {}

        out[key + '_disp'] = Struct(name='output_data', mode='vertex',
                                    data=dofs[:, :3], var_name=var_name,
                                    dofs=dof_names[:3])

        out[key + '_rot'] = Struct(name='output_data', mode='vertex',
                                   data=dofs[:, 3:], var_name=var_name,
                                   dofs=dof_names[3:])

        return out

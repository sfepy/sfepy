"""
Dirichlet, periodic and linear combination boundary condition classes.
"""
from sfepy.base.base import *

class BoundaryConditions(Container):
    """
    Container for boundary conditions.
    """
    @staticmethod
    def from_conf(conf):
        bcs = []
        for key, bc in conf.iteritems():
            if 'ebc' in key:
                bc = EssentialBC(bc.name, bc.region, bc.dofs, key=key)

            elif 'epbc' in key:
                bc = PeriodicBC(bc.name, bc.region, bc.dofs, bc.match, key=key)
    
            elif 'lcbc' in key:
                bc = LinearCombinationBC(bc.name, bc.region, bc.dofs, key=key)

            else:
                raise ValueError('unknown BC type! (%s)' % key)

            bcs.append(bc)

        obj = BoundaryConditions(bcs)
        return obj

    def group_by_variables(self, groups=None):
        """
        Group boundary conditions of each variable.
        """
        if groups is None:
            out = {}

        else:
            out = groups

        for bc in self:
            for single_bc in bc.iter_single():
                vname = single_bc.dofs[0].split( '.' )[0]
                out.setdefault(vname, BoundaryConditions()).append(single_bc)

        return out

    def canonize_dof_names(self, dofs):
        """
        Canonize the DOF names using the full list of DOFs of a
        variable.
        """
        for bc in self:
            bc.canonize_dof_names(dofs)

    def sort(self):
        """
        Sort boundary conditions by their key.
        """
        self._objs.sort(cmp=lambda i1, i2: cmp(i1.key, i2.key))
        self.update()

class BoundaryCondition(Struct):
    """
    Common boundary condition methods.
    """
    def iter_single(self):
        for dofs, val in self.dofs.iteritems():
            single_bc = self.copy()
            single_bc.dofs = [dofs, val]
            yield single_bc

    def canonize_dof_names(self, dofs):
        """
        Canonize the DOF names using the full list of DOFs of a
        variable.
        
        Assumes single BC instance.
        """
        vname, dd = self.dofs[0].split( '.' )
        if dd == 'all':
            cdofs = dofs
        elif dd[0] == '[':
            cdofs = [vname + '.' + ii.strip()
                     for ii in dd[1:-1].split( ',' )]
        else:
            cdofs = [self.dofs[0]]

        self.dofs[0] = cdofs

    def get_dofs(self):
        return self.cdofs

class EssentialBC(BoundaryCondition):
    """
    Essential boundary condidion.

    Parameters
    ----------
    name : str
        The boundary condition name.
    region : Region instance
        The region where the boundary condition is applied.
    dofs : dict
        The boundary condition specification defining the constrained
        DOFs and their values.
    key : str, optional
        The sorting key.
    """
    def __init__(self, name, region, dofs, key=''):
        Struct.__init__(self, name=name, region=region, dofs=dofs, key=key)

class PeriodicBC(BoundaryCondition):
    """
    Periodic boundary condidion.

    Parameters
    ----------
    name : str
        The boundary condition name.
    regions : list of two Region instances
        The master region and the slave region where the DOFs should match.
    dofs : dict
        The boundary condition specification defining the DOFs in the master
        region and the corresponding DOFs in the slave region.
    match : Function instance
        The function that should determine corresponding nodes in the
        two regions.
    key : str, optional
        The sorting key.
    """
    def __init__(self, name, regions, dofs, match, key=''):
        Struct.__init__(self, name=name, regions=regions, dofs=dofs,
                        match=match, key=key)


class LinearCombinationBC(BoundaryCondition):
    """
    Linear combination boundary condidion.

    Parameters
    ----------
    name : str
        The boundary condition name.
    region : Region instance
        The region where the boundary condition is applied.
    dofs : dict
        The boundary condition specification defining the constrained
        DOFs and the contraint type.
    key : str, optional
        The sorting key.
    """
    def __init__(self, name, region, dofs, key=''):
        Struct.__init__(self, name=name, region=region, dofs=dofs, key=key)

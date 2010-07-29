"""
The Dirichlet, periodic and linear combination boundary condition
classes, as well as the initial condition class.
"""
from sfepy.base.base import *

def _get_region(name, regions, bc_name):
    try:
        region = regions[name]
    except IndexError:
        msg = "no region '%s' used in condition %s!" % (name, bc_name)
        raise IndexError( msg )

    return region

class Conditions(Container):
    """
    Container for various conditions.
    """
    @staticmethod
    def from_conf(conf, regions):
        conds = []
        for key, cc in conf.iteritems():
            if 'ebc' in key:
                region = _get_region(cc.region, regions, cc.name)
                cond = EssentialBC(cc.name, region, cc.dofs, key=key)

            elif 'epbc' in key:
                rs = [_get_region(ii, regions, cc.name) for ii in cc.region]
                cond = PeriodicBC(cc.name, rs, cc.dofs, cc.match, key=key)

            elif 'lcbc' in key:
                region = _get_region(cc.region, regions, cc.name)
                cond = LinearCombinationBC(cc.name, region, cc.dofs, key=key)

            elif 'ic' in key:
                region = _get_region(cc.region, regions, cc.name)
                cond = InitialCondition(cc.name, region, cc.dofs, key=key)

            else:
                raise ValueError('unknown condition type! (%s)' % key)

            conds.append(cond)

        obj = Conditions(conds)
        return obj

    def group_by_variables(self, groups=None):
        """
        Group boundary conditions of each variable. Each condition is a
        group is a single condition.

        Parameters
        ----------
        groups : dict, optional
            If present, update the `groups` dictionary.

        Returns
        -------
        out : dict
            The dictionary with variable names as keys and lists of
            single condition instances as values.
        """
        if groups is None:
            out = {}

        else:
            out = groups

        for cond in self:
            for single_cond in cond.iter_single():
                vname = single_cond.dofs[0].split('.')[0]
                out.setdefault(vname, Conditions()).append(single_cond)

        return out

    def canonize_dof_names(self, dofs):
        """
        Canonize the DOF names using the full list of DOFs of a
        variable.
        """
        for cond in self:
            cond.canonize_dof_names(dofs)

    def sort(self):
        """
        Sort boundary conditions by their key.
        """
        self._objs.sort(cmp=lambda i1, i2: cmp(i1.key, i2.key))
        self.update()

def _canonize(dofs, all_dofs):
    """
    Helper function.
    """
    vname, dd = dofs.split('.')

    if dd == 'all':
        cdofs = all_dofs

    elif dd[0] == '[':
        cdofs = [vname + '.' + ii.strip()
                 for ii in dd[1:-1].split(',')]

    else:
        cdofs = [dofs]

    return cdofs

class Condition(Struct):
    """
    Common boundary condition methods.
    """
    def __init__(self, name, **kwargs):
        Struct.__init__(self, name=name, **kwargs)
        self.is_single = False

    def iter_single(self):
        """
        Create a single condition instance for each item in self.dofs
        and yield it.
        """
        for dofs, val in self.dofs.iteritems():
            single_cond = self.copy(name=self.name)
            single_cond.is_single = True
            single_cond.dofs = [dofs, val]
            yield single_cond

    def canonize_dof_names(self, dofs):
        """
        Canonize the DOF names using the full list of DOFs of a
        variable.
        
        Assumes single condition instance.
        """
        self.dofs[0] = _canonize(self.dofs[0], dofs)

class EssentialBC(Condition):
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
        Condition.__init__(self, name=name, region=region, dofs=dofs, key=key)

class PeriodicBC(Condition):
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
    match : str
        The name of function for matching corresponding nodes in the
        two regions.
    key : str, optional
        The sorting key.
    """
    def __init__(self, name, regions, dofs, match, key=''):
        Condition.__init__(self, name=name, regions=regions, dofs=dofs,
                           match=match, key=key)

    def canonize_dof_names(self, dofs):
        """
        Canonize the DOF names using the full list of DOFs of a
        variable.
        
        Assumes single condition instance.
        """
        self.dofs[0] = _canonize(self.dofs[0], dofs)
        self.dofs[1] = _canonize(self.dofs[1], dofs)

class LinearCombinationBC(Condition):
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
        DOFs and the constraint type.
    key : str, optional
        The sorting key.
    """
    def __init__(self, name, region, dofs, key=''):
        Condition.__init__(self, name=name, region=region, dofs=dofs, key=key)

class InitialCondition(Condition):
    """
    Initial condidion.

    Parameters
    ----------
    name : str
        The initial condition name.
    region : Region instance
        The region where the initial condition is applied.
    dofs : dict
        The initial condition specification defining the constrained
        DOFs and their values.
    key : str, optional
        The sorting key.
    """
    def __init__(self, name, region, dofs, key=''):
        Condition.__init__(self, name=name, region=region, dofs=dofs, key=key)

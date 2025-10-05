"""
The Dirichlet, periodic and linear combination boundary condition
classes, as well as the initial condition class.
"""
import numpy as nm

from sfepy.base.base import Container, Struct, is_sequence
from sfepy.discrete.functions import Function

def get_condition_value(val, functions, kind, name):
    """
    Check a boundary/initial condition value type and return the value or
    corresponding function.
    """
    if type(val) == str:
        if functions is not None:
            try:
                fun = functions[val]

            except IndexError:
                raise ValueError('unknown function %s given for %s %s!'
                                 % (val, kind, name))

        else:
            raise ValueError('no functions given for %s %s!' % (kind, name))

    elif (isinstance(val, Function) or nm.isscalar(val)
          or isinstance(val, nm.ndarray)):
        fun = val

    elif is_sequence(val):
        fun = nm.array(val)

    else:
        raise ValueError('unknown value type for %s %s!'
                         % (kind, name))

    return fun

def _get_region(name, regions, bc_name):
    try:
        region = regions[name]
    except IndexError:
        msg = "no region '%s' used in condition %s!" % (name, bc_name)
        raise IndexError(msg)

    return region

class Conditions(Container):
    """
    Container for various conditions.
    """
    @staticmethod
    def from_conf(conf, regions):
        conds = []
        for key, cc in conf.items():
            times = cc.get('times', None)



            if key.startswith("ebc"):
                region = _get_region(cc.region, regions, cc.name)
                cond = EssentialBC(cc.name, region, cc.dofs, key=key,
                                   times=times)

            elif key.startswith("epbc"):
                rs = [_get_region(ii, regions, cc.name) for ii in cc.region]
                cond = PeriodicBC(cc.name, rs, cc.dofs, cc.match, key=key,
                                   times=times)

            elif key.startswith('lcbc'):
                if isinstance(cc.region, str):
                    rs = [_get_region(cc.region, regions, cc.name), None]

                else:
                    rs = [_get_region(ii, regions, cc.name)
                          for ii in cc.region]

                cond = LinearCombinationBC(cc.name, rs, cc.dofs,
                                           cc.dof_map_fun, cc.kind,
                                           key=key,
                                           times=times,
                                           arguments=cc.get('arguments', None))

            elif key.startswith('dgebc'):
                region = _get_region(cc.region, regions, cc.name)
                cond = DGEssentialBC(cc.name, region, cc.dofs, key=key,
                                     times=times)

            elif key.startswith('dgepbc'):
                rs = [_get_region(ii, regions, cc.name) for ii in cc.region]
                cond = DGPeriodicBC(cc.name, rs, cc.dofs, cc.match, key=key,
                                    times=times)

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
        self._objs.sort(key=lambda a: a.key)
        self.update()

    def zero_dofs(self):
        """
        Set all boundary condition values to zero, if applicable.
        """
        for cond in self:
            if isinstance(cond, EssentialBC):
                cond.zero_dofs()

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
        for dofs, val in self.dofs.items():
            single_cond = self.copy(name=self.name)
            single_cond.is_single = True
            if 'grad' in dofs:
                # extract variable name from grad.<var-name>.all dofs
                dofs = ".".join((dofs.split(".")[1:]))
                # mark condition as diff
                single_cond.diff = 1

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
    times : list or str, optional
        The list of time intervals or a function returning True at time
        steps, when the condition applies.
    """
    def __init__(self, name, region, dofs, key='', times=None):
        Condition.__init__(self, name=name, region=region, dofs=dofs, key=key,
                           times=times)

    def zero_dofs(self):
        """
        Set all essential boundary condition values to zero.
        """
        if self.is_single:
            self.dofs[1] = 0.0

        else:
            new_dofs = {}
            for key in self.dofs.keys():
                new_dofs[key] = 0.0

            self.dofs = new_dofs

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
    times : list or str, optional
        The list of time intervals or a function returning True at time
        steps, when the condition applies.
    """
    def __init__(self, name, regions, dofs, match, key='', times=None):
        Condition.__init__(self, name=name, regions=regions, dofs=dofs,
                           match=match, key=key, times=times)

    def canonize_dof_names(self, dofs):
        """
        Canonize the DOF names using the full list of DOFs of a
        variable.

        Assumes single condition instance.
        """
        self.dofs[0] = _canonize(self.dofs[0], dofs)
        self.dofs[1] = _canonize(self.dofs[1], dofs)

class DGPeriodicBC(PeriodicBC):
    """
    This class is empty, it serves the same purpose
    as PeriodicBC, and is created only for branching in
    dof_info.py
    """
    pass

class DGEssentialBC(EssentialBC):
    """
    This class is empty, it serves the same purpose
    as EssentialBC, and is created only for branching in
    dof_info.py
    """

    def __init__(self, *args, diff=0, **kwargs):
        EssentialBC.__init__(self, *args, **kwargs)
        self.diff = diff

class LinearCombinationBC(Condition):
    """
    Linear combination boundary condidion.

    Parameters
    ----------
    name : str
        The boundary condition name.
    regions : list of two Region instances
        The constrained (master) DOFs region and the new (slave) DOFs
        region. The latter can be None if new DOFs are not field variable DOFs.
    dofs : dict
        The boundary condition specification defining the constrained
        DOFs and the new DOFs (can be None).
    dof_map_fun : str
        The name of function for mapping the constrained DOFs to new DOFs (can
        be None).
    kind : str
        The linear combination condition kind.
    key : str, optional
        The sorting key.
    times : list or str, optional
        The list of time intervals or a function returning True at time
        steps, when the condition applies.
    arguments: tuple, optional
        Additional arguments, depending on the condition kind.
    """
    def __init__(self, name, regions, dofs, dof_map_fun, kind, key='',
                 times=None, arguments=None):
        Condition.__init__(self, name=name, regions=regions, dofs=dofs,
                           dof_map_fun=dof_map_fun, kind=kind,
                           key=key, times=times, arguments=arguments)

    def get_var_names(self):
        """
        Get names of variables corresponding to the constrained and new DOFs.
        """
        names = [self.dofs[0].split('.')[0]]
        if self.dofs[1] is not None:
            names.append(self.dofs[1].split('.')[0])

        return names

    def canonize_dof_names(self, dofs0, dofs1=None):
        """
        Canonize the DOF names using the full list of DOFs of a
        variable.

        Assumes single condition instance.
        """
        self.dofs[0] = _canonize(self.dofs[0], dofs0)

        if self.dofs[1] is not None:
            self.dofs[1] = _canonize(self.dofs[1], dofs1)

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

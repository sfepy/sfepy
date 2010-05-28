"""
Classes holding information on global DOFs and mapping of all DOFs -
equations (active DOFs).
"""
from sfepy.base.base import *

class DofInfo(Struct):
    """
    Global DOF information, i.e. ordering of DOFs of the state (unknown)
    variables in the global state vector.
    """

    def __init__(self, name):
        Struct.__init__(self, name=name)

        self.n_var = 0
        self.var_names = []
        self.n_dof = {}
        self.ptr = [0]
        self.indx = {}
        self.details = {}

    def append_variable(self, var, active=False):
        """
        Append DOFs of the given variable.

        Parameters
        ----------
        var : Variable instance
            The variable to append.
        active : bool, optional
            When True, only active (non-constrained) DOFs are considered. 
        """
        name = var.name
        if name in self.var_names:
            raise ValueError('variable %s already present!' % name)

        self.var_names.append(name)

        self.n_dof[name], self.details[name] = var.get_dof_info(active=active)

        self.ptr.append(self.ptr[-1] + self.n_dof[name])

        ii = self.n_var
        self.indx[name] = slice(int(self.ptr[ii] ), int(self.ptr[ii+1]))
        
        self.n_var += 1

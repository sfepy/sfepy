"""
Module for handling state variables.
"""
from sfepy.base.base import Struct

class State(Struct):
    """
    Class holding/manipulating the state variables and corresponding DOF
    vectors.

    Manipulating the state class changes the underlying variables, and
    hence also the corresponding equations/terms (if any).

    Notes
    -----
    This class allows working with LCBC conditions in time-dependent
    problems, as it keeps track of the reduced DOF vector that cannot
    be reconstructed from the full DOF vector by using the usual
    `variables.strip_state_vector()`.
    """

    def __init__(self, variables, vec=None):
        """
        Create a State instance for the given variables.

        Parameters
        ----------
        variables : Variables instance
            The variables.
        vec : array, optional
            The (initial) DOF vector corresponding to the variables.
        """
        Struct.__init__(self, variables=variables, vec=vec, r_vec=None)

        if self.vec is None:
            self.vec = variables.create_state_vector()

        self.variables.set_data(self.vec)

    def copy(self, deep=False):
        """
        Copy the state. By default, the new state contains the same
        variables, and creates new DOF vectors. If `deep` is True, also
        the DOF vectors are copied.
        """
        if deep:
            other = State(self.variables, self.vec.copy())
            if self.r_vec is not None:
                other.r_vec = self.r_vec.copy()

        else:
            other = State(self.variables)

        return other

    def fill(self, value):
        """
        Fill the DOF vector with given value.
        """
        if self.r_vec is not None:
            self.r_vec.fill(value)

        self.vec.fill(value)

    def init_history(self):
        """
        Initialize variables with history.
        """
        self.variables.init_history()

    def apply_ebc(self, force_values=None):
        """
        Apply essential (Dirichlet) boundary conditions to the state.
        """
        self.variables.apply_ebc(self.vec, force_values=force_values)

    def has_ebc(self):
        """
        Test whether the essential (Dirichlet) boundary conditions have
        been applied to the DOF vector.
        """
        return self.variables.has_ebc(self.vec)

    def apply_ic(self, force_values=None):
        """
        Apply initial conditions to the state.
        """
        if self.r_vec is not None:
            raise ValueError('cannot re-apply initial conditions with LCBCs!')

        self.variables.apply_ic(self.vec, force_values=force_values)

    def get_reduced(self):
        """
        Get the reduced DOF vector, with EBC and PBC DOFs removed.
        """
        if self.variables.has_lcbc:
            if self.r_vec is None:
                r_vec = self.variables.strip_state_vector(self.vec)
                r_vec = self.variables.op_lcbc.T *  r_vec

            else:
                r_vec = self.r_vec

        else:
            r_vec = self.variables.strip_state_vector(self.vec)

        return r_vec

    def set_reduced(self, r_vec):
        """
        Set the reduced DOF vector, with EBC and PBC DOFs removed.
        """
        self.vec = self.variables.make_full_vec(r_vec)

        if self.variables.has_lcbc:
            self.r_vec = r_vec

    def set_full(self, vec, var_name=None):
        """
        Set the full DOF vector (including EBC and PBC DOFs). If
        `var_name` is given, set only the DOF sub-vector corresponding to
        the given variable.
        """
        if var_name is None:
            if self.variables.has_lcbc:
                raise ValueError('cannot set full DOF vector with LCBCs!')

            self.vec = vec

        else:
            var = self.variables[var_name]

            if var.has_lcbc:
                raise ValueError('cannot set full DOF vector with LCBCs!')

            self.variables.set_state_part(self.vec, vec, var_name)

    def __call__(self, var_name=None):
        """
        Get the full DOF vector (including EBC and PBC DOFs). If
        `var_name` is given, return only the DOF vector corresponding to
        the given variable.
        """
        if var_name is None:
            out = self.vec

        else:
            out = self.variables.get_state_part_view(self.vec, var_name)

        return out

    def get_parts(self):
        """
        Return parts of the DOF vector corresponding to individual state
        variables.

        Returns
        -------
        out : dict
            The dictionary of the DOF vector parts.
        """
        return self.variables.get_state_parts(self.vec)

    def create_output_dict(self, fill_value=None, var_info=None,
                           extend=True):
        """
        Transforms state to an output dictionary, that can be
        passed as 'out' kwarg to Mesh.write().

        Then the dictionary entries are formed by components of the
        state vector corresponding to the unknown variables, each
        transformed to shape (n_mesh_nod, n_dof per node) - all values
        in extra (higher order) nodes are removed.

        Examples
        --------
        >>> out = state.create_output_dict()
        >>> problem.save_state('file.vtk', out=out)
        """
        return self.variables.state_to_output(self.vec, fill_value,
                                              var_info, extend)

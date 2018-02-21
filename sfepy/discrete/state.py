"""
Module for handling state variables.
"""
from __future__ import absolute_import
import numpy as nm

from sfepy.base.base import Struct
import six

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

    @staticmethod
    def from_variables(variables):
        """
        Create a State instance for the given variables.

        The DOF vector is created using the DOF data in `variables`.

        Parameters
        ----------
        variables : Variables instance
            The variables.
        """
        parts = variables.get_state_parts()
        vec = variables.create_state_vector()

        for key, part in six.iteritems(parts):
            indx = variables.get_indx(key)
            vec[indx] = part

        return State(variables, vec)

    def __init__(self, variables, vec=None, preserve_caches=False):
        """
        Create a State instance for the given variables.

        Parameters
        ----------
        variables : Variables instance
            The variables.
        vec : array, optional
            The (initial) DOF vector corresponding to the variables.
        preserve_caches : bool
            If True, do not invalidate evaluate caches of variables.
        """
        Struct.__init__(self, variables=variables, vec=vec, r_vec=None)

        if self.vec is None:
            self.vec = variables.create_state_vector()

        self.variables.set_data(self.vec, preserve_caches=preserve_caches)

    def copy(self, deep=False, preserve_caches=False):
        """
        Copy the state. By default, the new state contains the same
        variables, and creates new DOF vectors. If `deep` is True, also
        the DOF vectors are copied.

        Parameters
        ----------
        deep : bool
            If True, make a copy of the DOF vectors.
        preserve_caches : bool
            If True, do not invalidate evaluate caches of variables.
        """
        if deep:
            other = State(self.variables, self.vec.copy(), preserve_caches=True)
            if self.r_vec is not None:
                other.r_vec = self.r_vec.copy()

        else:
            other = State(self.variables, preserve_caches=True)

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

    def get_reduced(self, follow_epbc=False):
        """
        Get the reduced DOF vector, with EBC and PBC DOFs removed.
        """
        strip = self.variables.strip_state_vector

        if self.variables.has_lcbc:
            if self.r_vec is None:
                r_vec = strip(self.vec, follow_epbc=follow_epbc)
                # This just sets the correct vector size (wrong values)!
                r_vec = self.variables.mtx_lcbc.T * r_vec

            else:
                r_vec = self.r_vec

        else:
            r_vec = strip(self.vec, follow_epbc=follow_epbc)

        return r_vec

    def set_reduced(self, r_vec, preserve_caches=False):
        """
        Set the reduced DOF vector, with EBC and PBC DOFs removed.

        Parameters
        ----------
        r_vec : array
            The reduced DOF vector corresponding to the variables.
        preserve_caches : bool
            If True, do not invalidate evaluate caches of variables.
        """
        self.vec = self.variables.make_full_vec(r_vec)

        if self.variables.has_lcbc:
            self.r_vec = r_vec

        self.variables.set_data(self.vec, preserve_caches=preserve_caches)

    def set_full(self, vec, var_name=None, force=False):
        """
        Set the full DOF vector (including EBC and PBC DOFs). If
        `var_name` is given, set only the DOF sub-vector corresponding
        to the given variable. If `force` is True, setting variables
        with LCBC DOFs is allowed.
        """
        if var_name is None:
            if self.variables.has_lcbc and not force:
                raise ValueError('cannot set full DOF vector with LCBCs!')

            self.vec = vec
            self.variables.set_data(self.vec)

        else:
            var = self.variables[var_name]

            if var.has_lcbc and not force:
                raise ValueError('cannot set full DOF vector with LCBCs!')

            self.variables.set_state_part(self.vec, vec, var_name)
            var.set_data(self.vec, self.variables.get_indx(var_name))

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

    def set_parts(self, parts, force=False):
        """
        Set parts of the DOF vector corresponding to individual state
        variables.

        Parameters
        ----------
        parts : dict
            The dictionary of the DOF vector parts.
        """
        if self.variables.has_lcbc and not force:
            raise ValueError('cannot set full DOF vector with LCBCs!')

        self.variables.set_data(parts)
        for key, part in six.iteritems(parts):
            indx = self.variables.get_indx(key)
            self.vec[indx] = part

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

    def get_vec(self, active_only):
        if active_only:
            vec = self.get_reduced()

        else:
            vec = self()

        return vec

    def set_vec(self, vec, active_only):
        if active_only:
            self.set_reduced(vec, preserve_caches=True)

        else:
            self.set_full(vec)

    def create_output_dict(self, fill_value=None, var_info=None,
                           extend=True, linearization=None):
        """
        Transforms state to an output dictionary, that can be
        passed as 'out' kwarg to Mesh.write().

        Then the dictionary entries are formed by components of the
        state vector corresponding to unknown variables according to
        kind of linearization given by `linearization`.

        Examples
        --------
        >>> out = state.create_output_dict()
        >>> problem.save_state('file.vtk', out=out)
        """
        return self.variables.state_to_output(self.vec, fill_value,
                                              var_info, extend,
                                              linearization=linearization)

    def get_weighted_norm(self, vec, weights=None, return_weights=False):
        """
        Return the weighted norm of DOF vector `vec`.

        By default, each component of `vec` is weighted by the 1/norm of the
        corresponding state part, or 1 if the norm is zero. Alternatively, the
        weights can be provided explicitly using `weights` argument.

        Parameters
        ----------
        vec : array
            The DOF vector corresponding to the variables.
        weights : dict, optional
            If given, the weights are used instead of the norms of the state
            parts. Keys of the dictionary must be equal to the names of
            variables comprising the DOF vector.
        return_weights: bool
            If True, return also the used weights.

        Returns
        -------
        norm : float
            The weighted norm.
        weights : dict, optional
            If `return_weights` is True, the used weights.

        Examples
        --------
        >>> err = state0.get_weighted_norm(state() - state0())
        """
        if weights is None:
            parts = self.get_parts()

            weights = {}
            for key, part in six.iteritems(parts):
                pnorm = nm.linalg.norm(part)
                if pnorm < 10.0 * nm.finfo(nm.float64).eps:
                    pnorm = 1.0
                weights[key] = 1.0 / pnorm

        else:
            if set(weights.keys()) != self.variables.state:
                raise ValueError('weights keys have to be in %s!'
                                 % self.variables.state)

        wvec = vec.copy()
        for key in six.iterkeys(weights):
            indx = self.variables.get_indx(key)
            wvec[indx] *= weights[key]

        norm = nm.linalg.norm(wvec)

        if return_weights:
            return norm, weights

        else:
            return norm

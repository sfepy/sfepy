"""
Classes holding information on global DOFs and mapping of all DOFs -
equations (active DOFs).

Helper functions for the equation mapping.
"""
import numpy as nm
import scipy.sparse as sp

from sfepy.base.base import assert_, Struct, basestr
from sfepy.discrete.functions import Function
from sfepy.discrete.conditions import get_condition_value, EssentialBC, \
    PeriodicBC, DGPeriodicBC, DGEssentialBC


def expand_nodes_to_dofs(nods, n_dof_per_node):
    """
    Expand DOF node indices into DOFs given a constant number of DOFs
    per node.
    """
    dofs = nm.repeat(nods, n_dof_per_node)
    dofs.shape = (nods.shape[0], n_dof_per_node)

    idof = nm.arange(n_dof_per_node, dtype=nm.int32)

    dofs = n_dof_per_node * dofs + idof

    return dofs

def expand_nodes_to_equations(nods, dof_names, all_dof_names):
    """
    Expand vector of node indices to equations (DOF indices) based on
    the DOF-per-node count.

    DOF names must be already canonized.

    Returns
    -------
    eq : array
        The equations/DOF indices in the node-by-node order.
    """
    dpn = len(all_dof_names)
    nc = len(dof_names)

    eq = nm.empty(len(nods) * nc, dtype=nm.int32)
    for ii, dof in enumerate(dof_names):
        idof = all_dof_names.index(dof)
        eq[ii::nc] = dpn * nods + idof
    return eq

def resolve_chains(master_slave, chains):
    """
    Resolve EPBC chains - e.g. in corner nodes.
    """
    for chain in chains:
        slave = chain[-1]
        master_slave[chain[:-1]] = slave + 1
        master_slave[slave] = - chain[0] - 1 # Any of masters...

def group_chains(chain_list):
    """
    Group EPBC chains.
    """
    chains = []
    while len(chain_list):
        chain = set(chain_list.pop(0))
        ## print ':', chain
        ii = 0
        while ii < len(chain_list):
            c1 = sorted(chain_list[ii])
            ## print '--', ii, c1, chain
            is0 = c1[0] in chain
            is1 = c1[1] in chain

            if is0 and is1:
                chain_list.pop(ii)
            elif is0 or is1:
                chain.update(c1)
                chain_list.pop(ii)
                ii = 0
            else:
                ii += 1
            ## print ii, chain, chain_list
        ## print '->', chain
        ## print chain_list

        chains.append(list(chain))

    ## print 'EPBC chain groups:', chains
    aux = {}
    for chain in chains:
        aux.setdefault(len(chain), [0])[0] += 1
    ## print 'EPBC chain counts:', aux

    return chains

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

    def _update_after_append(self, name):
        self.ptr.append(self.ptr[-1] + self.n_dof[name])

        ii = self.n_var
        self.indx[name] = slice(int(self.ptr[ii]), int(self.ptr[ii+1]))

        self.n_var += 1

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
        self._update_after_append(name)

    def append_raw(self, name, n_dof):
        """
        Append raw DOFs.

        Parameters
        ----------
        name : str
            The name of variable the DOFs correspond to.
        n_dof : int
            The number of DOFs.
        """
        if name in self.var_names:
            raise ValueError('variable %s already present!' % name)

        self.var_names.append(name)

        self.n_dof[name], self.details[name] = n_dof, None
        self._update_after_append(name)

    def update(self, name, n_dof):
        """
        Set the number of DOFs of the given variable.

        Parameters
        ----------
        name : str
            The name of variable the DOFs correspond to.
        n_dof : int
            The number of DOFs.
        """
        if not name in self.var_names:
            raise ValueError('variable %s is not present!' % name)

        ii = self.var_names.index(name)
        delta = n_dof - self.n_dof[name]

        self.n_dof[name] = n_dof

        for iv, nn in enumerate(self.var_names[ii:]):
            self.ptr[ii+iv+1] += delta
            self.indx[nn] = slice(self.ptr[ii+iv], self.ptr[ii+iv+1])

    def get_info(self, var_name):
        """
        Return information on DOFs of the given variable.

        Parameters
        ----------
        var_name : str
            The name of the variable.
        """
        return Struct(name='%s_dof_info' % var_name,
                      var_name=var_name,
                      n_dof=self.n_dof[var_name],
                      indx=self.indx[var_name],
                      details=self.details[var_name])

    def get_subset_info(self, var_names):
        """
        Return global DOF information for selected variables
        only. Silently ignores non-existing variable names.

        Parameters
        ----------
        var_names : list
            The names of the selected variables.
        """
        di = DofInfo(self.name + ':subset')
        for var_name in var_names:
            if var_name not in self.var_names:
                continue

            di.append_raw(var_name, self.n_dof[var_name])

        return di

    def get_n_dof_total(self):
        """
        Return the total number of DOFs of all state variables.
        """
        return self.ptr[-1]

def is_active_bc(bc, ts=None, functions=None):
    """
    Check whether the given boundary condition is active in the current
    time.

    Returns
    -------
    active : bool
        True if the condition `bc` is active.
    """
    if (bc.times is None) or (ts is None):
        active = True

    elif isinstance(bc.times, list):
        for tt in bc.times:
            if tt[0] <= ts.time < tt[1]:
                active = True
                break

        else:
            active = False

    else:
        if isinstance(bc.times, basestr):
            if functions is not None:
                fun = functions[bc.times]

            else:
                raise ValueError('no functions given for bc %s!' % bc.name)

        elif isinstance(bc.times, Function):
            fun = bc.times

        else:
            raise ValueError('unknown times type! (%s)'
                             % type(bc.times))

        active = fun(ts)

    return active

class EquationMap(Struct):
    """
    Map all DOFs to equations for active DOFs.
    """

    def __init__(self, name, dof_names, var_di):
        Struct.__init__(self, name=name, dof_names=dof_names, var_di=var_di)

        self.dpn = len(self.dof_names)
        self.eq = nm.arange(var_di.n_dof, dtype=nm.int32)

        self.n_dg_ebc = 0
        self.dg_ebc_names = {}
        self.dg_ebc = {}
        self.dg_ebc_val = {}

        self.n_dg_epbc = 0
        self.dg_epbc_names = []
        self.dg_epbc = []

    def _init_empty(self, field):
        self.val_ebc = nm.empty((0,), dtype=field.dtype)

        if field.get('unused_dofs') is None:
            self.eqi = nm.arange(self.var_di.n_dof, dtype=nm.int32)

        else:
            self._mark_unused(field)
            self.eqi = nm.compress(self.eq >= 0, self.eq)
            self.eq[self.eqi] = nm.arange(self.eqi.shape[0], dtype=nm.int32)

        self.eq_ebc = nm.empty((0,), dtype=nm.int32)

        self.master = nm.empty((0,), dtype=nm.int32)
        self.slave = nm.empty((0,), dtype=nm.int32)

        self.n_eq = self.eqi.shape[0]
        self.n_ebc = self.eq_ebc.shape[0]
        self.n_epbc = self.master.shape[0]

    def _mark_unused(self, field):
        unused_dofs = field.get('unused_dofs')
        if unused_dofs is not None:
            unused = expand_nodes_to_equations(field.unused_dofs,
                                               self.dof_names, self.dof_names)
            self.eq[unused] = -3

    def map_equations(self, bcs, field, ts, functions, problem=None,
                      warn=False):
        """
        Create the mapping of active DOFs from/to all DOFs.

        Parameters
        ----------
        bcs : Conditions instance
            The Dirichlet or periodic boundary conditions (single
            condition instances). The dof names in the conditions must
            already be canonized.
        field : Field instance
            The field of the variable holding the DOFs.
        ts : TimeStepper instance
            The time stepper.
        functions : Functions instance
            The registered functions.
        problem : Problem instance, optional
            The problem that can be passed to user functions as a context.
        warn : bool, optional
            If True, warn about BC on non-existent nodes.

        Returns
        -------
        active_bcs : set
            The set of boundary conditions active in the current time.

        Notes
        -----
        - Periodic bc: master and slave DOFs must belong to the same
          field (variables can differ, though).
        """
        if bcs is None:
            self._init_empty(field)
            return set()

        eq_ebc = nm.zeros((self.var_di.n_dof,), dtype=nm.int32)
        val_ebc = nm.zeros((self.var_di.n_dof,), dtype=field.dtype)
        master_slave = nm.zeros((self.var_di.n_dof,), dtype=nm.int32)
        chains = []

        active_bcs = set()
        for bc in bcs:
            # Skip conditions that are not active in the current time.
            if not is_active_bc(bc, ts=ts, functions=functions):
                continue

            active_bcs.add(bc.key)
            if isinstance(bc, DGEssentialBC):
                ntype = "DGEBC"
                region = bc.region
            elif isinstance(bc, DGPeriodicBC):
                ntype = "DGEPBC"
                region = bc.regions[0]
            elif isinstance(bc, EssentialBC):
                ntype = 'EBC'
                region = bc.region
            elif isinstance(bc, PeriodicBC):
                ntype = 'EPBC'
                region = bc.regions[0]

            if warn:
                clean_msg = ('warning: ignoring nonexistent %s node (%s) in '
                             % (ntype, self.var_di.var_name))
            else:
                clean_msg = None

            # Get master region nodes.
            master_nod_list = field.get_dofs_in_region(region)
            if len(master_nod_list) == 0:
                continue

            if ntype == 'EBC': # EBC.
                dofs, val = bc.dofs
                ##
                # Evaluate EBC values.
                fun = get_condition_value(val, functions, 'EBC', bc.name)
                if isinstance(fun, Function):
                    aux = fun
                    fun = lambda coors: aux(ts, coors,
                                            bc=bc, problem=problem)

                nods, vv = field.set_dofs(fun, region, len(dofs), clean_msg)

                eq = expand_nodes_to_equations(nods, dofs, self.dof_names)
                # Duplicates removed here...
                eq_ebc[eq] = 1
                if vv is not None: val_ebc[eq] = nm.ravel(vv)
            elif ntype == "DGEBC":

                dofs, val = bc.dofs
                ##
                # Evaluate EBC values.
                fun = get_condition_value(val, functions, 'EBC', bc.name)
                if isinstance(fun, Function):
                    aux = fun
                    fun = lambda coors: aux(ts, coors,
                                            bc=bc, problem=problem)

                values = field.get_bc_facet_values(fun, region, diff=bc.diff)
                bc2bfi = field.get_bc_facet_idx(region)

                self.dg_ebc_val.setdefault(bc.diff, []).append(values)
                self.dg_ebc.setdefault(bc.diff, []).append(bc2bfi)
                self.n_dg_ebc += 1
            elif ntype == "DGEPBC":

                # ensure matching boundaries?
                master_bc2bfi = field.get_bc_facet_idx(region)
                slave_bc2bfi = field.get_bc_facet_idx(bc.regions[1])

                self.dg_epbc.append((master_bc2bfi, slave_bc2bfi))
                self.n_dg_epbc += 1

            else: # EPBC.
                region = bc.regions[1]
                slave_nod_list = field.get_dofs_in_region(region)

                nmaster = nm.unique(nm.hstack(master_nod_list))
                # Treat fields not covering the whole domain.
                if nmaster[0] == -1:
                    nmaster = nmaster[1:]

                nslave = nm.unique(nm.hstack(slave_nod_list))
                # Treat fields not covering the whole domain.
                if nslave[0] == -1:
                    nslave = nslave[1:]

                ## print nmaster + 1
                ## print nslave + 1
                if nmaster.shape != nslave.shape:
                    msg = 'EPBC list lengths do not match!\n(%s,\n %s)' %\
                          (nmaster, nslave)
                    raise ValueError(msg)

                if (nmaster.shape[0] == 0) and (nslave.shape[0] == 0):
                    continue

                mcoor = field.get_coor(nmaster)
                scoor = field.get_coor(nslave)

                fun = get_condition_value(bc.match, functions, 'EPBC', bc.name)
                if isinstance(fun, Function):
                    i1, i2 = fun(mcoor, scoor)

                else:
                    i1, i2 = fun

                ## print nm.c_[mcoor[i1], scoor[i2]]
                ## print nm.c_[nmaster[i1], nslave[i2]] + 1

                meq = expand_nodes_to_equations(nmaster[i1], bc.dofs[0],
                                                self.dof_names)
                seq = expand_nodes_to_equations(nslave[i2], bc.dofs[1],
                                                self.dof_names)

                m_assigned = nm.where(master_slave[meq] != 0)[0]
                s_assigned = nm.where(master_slave[seq] != 0)[0]
                if m_assigned.size or s_assigned.size: # Chain EPBC.
                    aux = master_slave[meq[m_assigned]]
                    sgn = nm.sign(aux)
                    om_chain = zip(meq[m_assigned], (aux - sgn) * sgn)
                    chains.extend(om_chain)

                    aux = master_slave[seq[s_assigned]]
                    sgn = nm.sign(aux)
                    os_chain = zip(seq[s_assigned], (aux - sgn) * sgn)
                    chains.extend(os_chain)

                    m_chain = zip(meq[m_assigned], seq[m_assigned])
                    chains.extend(m_chain)

                    msd = nm.setdiff1d(s_assigned, m_assigned)
                    s_chain = zip(meq[msd], seq[msd])
                    chains.extend(s_chain)

                    msa = nm.union1d(m_assigned, s_assigned)
                    ii = nm.setdiff1d(nm.arange(meq.size), msa)
                    master_slave[meq[ii]] = seq[ii] + 1
                    master_slave[seq[ii]] = - meq[ii] - 1

                else:
                    master_slave[meq] = seq + 1
                    master_slave[seq] = - meq - 1

        chains = group_chains(chains)
        resolve_chains(master_slave, chains)

        ii = nm.argwhere(eq_ebc == 1)
        self.eq_ebc = nm.atleast_1d(ii.squeeze())
        self.val_ebc = nm.atleast_1d(val_ebc[ii].squeeze())
        # add axis in case we squeezed too hard
        self.master = nm.atleast_1d(nm.argwhere(master_slave > 0).squeeze())
        self.slave = master_slave[self.master] - 1

        assert_((self.eq_ebc.shape == self.val_ebc.shape))
        self.eq[self.eq_ebc] = -2
        self.eq[self.master] = -1

        self._mark_unused(field)

        self.eqi = nm.compress(self.eq >= 0, self.eq)
        self.eq[self.eqi] = nm.arange(self.eqi.shape[0], dtype=nm.int32)
        self.eq[self.master] = self.eq[self.slave]
        self.n_eq = self.eqi.shape[0]
        self.n_ebc = self.eq_ebc.shape[0]
        self.n_epbc = self.master.shape[0]

        return active_bcs

    def get_operator(self):
        """
        Get the matrix operator :math:`R` corresponding to the equation
        mapping, such that the restricted matrix :math:`A_r` can be
        obtained from the full matrix :math:`A` by :math:`A_r = R^T A
        R`. All the matrices are w.r.t. a single variables that uses
        this mapping.

        Returns
        -------
        mtx : coo_matrix
            The matrix :math:`R`.
        """
        # EBC.
        rows = self.eqi
        cols = nm.arange(self.n_eq, dtype=nm.int32)

        # EPBC.
        ic = self.eq[self.slave]
        ii = ic >= 0
        rows = nm.r_[rows, self.master[ii]]
        cols = nm.r_[cols, ic[ii]]

        ones = nm.ones(rows.shape[0], dtype=nm.float64)
        mtx = sp.coo_matrix((ones, (rows, cols)),
                            shape=(self.eq.shape[0], self.n_eq))

        return mtx

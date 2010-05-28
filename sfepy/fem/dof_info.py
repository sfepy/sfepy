"""
Classes holding information on global DOFs and mapping of all DOFs -
equations (active DOFs).

Helper functions for the equation mapping.
"""
from sfepy.base.base import *

def expand_nodes_to_equations(nods, dof_names, all_dof_names):
    """
    Expand vector of node indices to equations (DOF indices) based on
    the DOF-per-node count.

    DOF names must be already canonized.
    """
    dpn = len(all_dof_names)
    
    eq = nm.array([], dtype=nm.int32)
    for dof in dof_names:
        idof = all_dof_names.index(dof)
        eq = nm.concatenate((eq, dpn * nods + idof))
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
        ## pause()

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

    def get_info(self, var_name):
        """
        Return information on DOFs of the given variable.

        Parameters
        ----------
        var_name : str
            The name of the variable.
        """
        return Struct(name = '%s_dof_info' % var_name,
                      var_name = var_name,
                      n_dof = self.n_dof[var_name],
                      indx = self.indx[var_name],
                      details = self.details[var_name])

class EquationMap(Struct):
    """
    Map all DOFs to equations for active DOFs.
    """

    def __init__(self, name, dof_names, var_di):
        Struct.__init__(self, name=name, dof_names=dof_names, var_di=var_di)

        self.dpn = len(self.dof_names)
        self.eq = nm.arange(var_di.n_dof, dtype=nm.int32)

    def _init_empty(self):
        self.eqi = nm.arange(self.var_di.n_dof, dtype=nm.int32)
        self.eq_ebc = nm.empty((0,), dtype=nm.int32)

        self.n_eq = self.eqi.shape[0]
        self.n_ebc = self.eq_ebc.shape[0]

        self.master = nm.empty((0,), dtype=nm.int32)
        self.slave = nm.empty((0,), dtype=nm.int32)
        
    def map_equations(self, bcs, field, regions, ts, functions, warn=False):
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
        regions : list
            List of regions.
        ts : TimeStepper instance
            The time stepper.
        functions : Functions instance
            The registered functions.
        warn : bool, optional
            If True, warn about BC on non-existent nodes.

        Notes
        -----
        - Periodic bc: master and slave DOFs must belong to the same
          field (variables can differ, though).
        """
        if bcs is None:
            self.val_ebc = nm.empty((0,), dtype=field.dtype)
            self._init_empty()
            return

        eq_ebc = nm.zeros((self.var_di.n_dof,), dtype=nm.int32)
        val_ebc = nm.zeros((self.var_di.n_dof,), dtype=field.dtype)
        master_slave = nm.zeros((self.var_di.n_dof,), dtype=nm.int32)
        chains = []

        for bc in bcs:
            if 'ebc' in bc.key:
                ntype = 'EBC'
                rname = bc.region
            else:
                ntype = 'EPBC'
                rname = bc.regions[0]

            try:
                region = regions[rname]
            except IndexError:
                msg = "no region '%s' used in BC %s!" % (rname, bc)
                raise IndexError( msg )

            ## print ir, key, bc
            ## debug()

            if warn:
                clean_msg = ('warning: ignoring nonexistent' \
                             ' %s node (%s) in ' % (ntype, self.var_di.var_name))
            else:
                clean_msg = None

            # Get master region nodes.
            master_nod_list = region.get_field_nodes(field, clean=True,
                                                     warn=clean_msg)
            if len(master_nod_list) == 0:
                continue

            if ntype == 'EBC': # EBC.
                dofs, val = bc.dofs
                ##
                # Evaluate EBC values.
                nods = nm.unique(nm.hstack(master_nod_list))
                coor = field.get_coor(nods)

                if type(val) == str:
                    fun = functions[val]
                    vv = fun(ts, coor, bc=bc)

                else:
                    vv = nm.repeat([val], nods.shape[0] * len(dofs))

                eq = expand_nodes_to_equations(nods, dofs, self.dof_names)
                # Duplicates removed here...
                eq_ebc[eq] = 1
                if vv is not None: val_ebc[eq] = vv

            else: # EPBC.
                region = regions[bc.regions[1]]
                slave_nod_list = region.get_field_nodes(field, clean=True,
                                                        warn=clean_msg)
                ## print master_nod_list
                ## print slave_nod_list

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

                fun = functions[bc.match]
                i1, i2 = fun(mcoor, scoor)
               ## print nm.c_[mcoor[i1], scoor[i2]]
               ## print nm.c_[nmaster[i1], nslave[i2]] + 1

                meq = expand_nodes_to_equations(nmaster[i1], bc.dofs[0],
                                                self.dof_names)
                seq = expand_nodes_to_equations(nslave[i2], bc.dofs[1],
                                                self.dof_names)

                m_assigned = nm.where(master_slave[meq] != 0)[0]
                s_assigned = nm.where(master_slave[seq] != 0)[0]
                if m_assigned.size or s_assigned.size: # Chain EPBC.
                    ## print m_assigned, meq[m_assigned]
                    ## print s_assigned, seq[s_assigned]

                    aux = master_slave[meq[m_assigned]]
                    sgn = nm.sign(aux)
                    om_chain = zip(meq[m_assigned], (aux - sgn) * sgn)
                    ## print om_chain
                    chains.extend(om_chain)

                    aux = master_slave[seq[s_assigned]]
                    sgn = nm.sign(aux)
                    os_chain = zip(seq[s_assigned], (aux - sgn) * sgn)
                    ## print os_chain
                    chains.extend(os_chain)

                    m_chain = zip(meq[m_assigned], seq[m_assigned])
                    ## print m_chain
                    chains.extend(m_chain)

                    msd = nm.setdiff1d(s_assigned, m_assigned)
                    s_chain = zip(meq[msd], seq[msd])
                    ## print s_chain
                    chains.extend(s_chain)

                    msa = nm.union1d(m_assigned, s_assigned)
                    ii = nm.setdiff1d(nm.arange(meq.size), msa)
                    master_slave[meq[ii]] = seq[ii] + 1
                    master_slave[seq[ii]] = - meq[ii] - 1

                else:
                    master_slave[meq] = seq + 1
                    master_slave[seq] = - meq - 1
                ## print 'ms', master_slave
                ## print chains

        ## print master_slave
        chains = group_chains(chains)
        resolve_chains(master_slave, chains)

        ii = nm.argwhere(eq_ebc == 1)
        self.eq_ebc = nm.atleast_1d(ii.squeeze())
        self.val_ebc = nm.atleast_1d(val_ebc[ii].squeeze())
        self.master = nm.argwhere(master_slave > 0).squeeze()
        self.slave = master_slave[self.master] - 1

        assert_((self.eq_ebc.shape == self.val_ebc.shape))
        ## print self.eq_ebc.shape
        ## pause()
        self.eq[self.eq_ebc] = -2
        self.eq[self.master] = -1
        self.eqi = nm.compress(self.eq >= 0, self.eq)
        self.eq[self.eqi] = nm.arange(self.eqi.shape[0], dtype=nm.int32)
        self.eq[self.master] = self.eq[self.slave]
        self.n_eq = self.eqi.shape[0]
        self.n_ebc = self.eq_ebc.shape[0]
        self.n_epbc = self.master.shape[0]

"""
Classes holding information on global DOFs and mapping of all DOFs -
equations (active DOFs).

Helper functions for the equation mapping.
"""
from sfepy.base.base import *
from sfepy.fem.utils import compute_nodal_normals
from sfepy.fem.functions import Function
from sfepy.fem.conditions import EssentialBC

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
        
    def map_equations(self, bcs, field, ts, functions, warn=False):
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
            if isinstance(bc, EssentialBC):
                ntype = 'EBC'
                region = bc.region

            else:
                ntype = 'EPBC'
                region = bc.regions[0]

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

                elif isinstance(val, Function):
                    vv = val(ts, coor, bc=bc)

                else:
                    vv = nm.repeat([val], nods.shape[0] * len(dofs))

                eq = expand_nodes_to_equations(nods, dofs, self.dof_names)
                # Duplicates removed here...
                eq_ebc[eq] = 1
                if vv is not None: val_ebc[eq] = vv

            else: # EPBC.
                region = bc.regions[1]
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

class LCBCOperator(Struct):
    """
    Base class for LCBC operators.
    """

    def treat_pbcs(self, master_equations):
        """
        Treat dofs with periodic BC.
        """
        umeq, indx = nm.unique1d(master_equations, return_index=True)
        indx.sort()
        self.mtx = self.mtx[indx]

class RigidOperator(LCBCOperator):
    """
    Transformation matrix operator for rigid LCBCs.
    """

    def __init__(self, name, nodes, field, dof_names, all_dof_names):
        Struct.__init__(self, name=name, nodes=nodes, dof_names=dof_names)

        coors = field.get_coor(nodes)
        n_nod, dim = coors.shape

        mtx_e = nm.tile(nm.eye(dim, dtype=nm.float64), (n_nod, 1))

        if dim == 2:
            mtx_r = nm.empty((dim * n_nod, 1), dtype=nm.float64)
            mtx_r[0::dim,0] = -coors[:,1]
            mtx_r[1::dim,0] = coors[:,0]
            n_rigid_dof = 3

        elif dim == 3:
            mtx_r = nm.zeros((dim * n_nod, dim), dtype=nm.float64)
            mtx_r[0::dim,1] = coors[:,2]
            mtx_r[0::dim,2] = -coors[:,1]
            mtx_r[1::dim,0] = -coors[:,2]
            mtx_r[1::dim,2] = coors[:,0]
            mtx_r[2::dim,0] = coors[:,1]
            mtx_r[2::dim,1] = -coors[:,0]
            n_rigid_dof = 6

        else:
            msg = 'dimension in [2, 3]: %d' % dim
            raise ValueError(msg)

        self.n_dof = n_rigid_dof
        self.mtx = nm.hstack((mtx_r, mtx_e))

        # Strip unconstrained dofs.
        aux = dim * nm.arange(n_nod)
        indx = [aux + all_dof_names.index(dof) for dof in dof_names]
        indx = nm.array(indx).T.ravel()

        self.mtx = self.mtx[indx]

def _save_normals(filename, normals, region, mesh):
        nn = nm.zeros_like(mesh.coors)
        nmax = region.all_vertices.shape[0]
        nn[region.all_vertices] = normals[:nmax]
        out = {'normals' : Struct(name = 'output_data',
                                  mode = 'vertex', data = nn)}
        mesh.write(filename, out=out, io='auto')

class NoPenetrationOperator(LCBCOperator):
    """
    Transformation matrix operator for no-penetration LCBCs.
    """
    def __init__(self, name, nodes, region, field, dof_names, filename=None):
        Struct.__init__(self, name=name, nodes=nodes, dof_names=dof_names)

        dim = field.shape[0]
        assert_(len(dof_names) == dim)

        normals = compute_nodal_normals(nodes, region, field)

        if filename is not None:
            _save_normals(filename, normals, region, field.domain.mesh)

        ii = nm.abs(normals).argmax(1)
        n_nod, dim = normals.shape

        irs = set(range(dim))

        data = []
        rows = []
        cols = []
        for idim in xrange(dim):
            ic = nm.where(ii == idim)[0]
            if len(ic) == 0: continue
            ## print ic
            ## print idim

            ir = list(irs.difference([idim]))
            nn = nm.empty((len(ic), dim - 1), dtype=nm.float64)
            for ik, il in enumerate(ir):
                nn[:,ik] = - normals[ic,il] / normals[ic,idim]

            irn = dim * ic + idim
            ics = [(dim - 1) * ic + ik for ik in xrange(dim - 1)]
            for ik in xrange(dim - 1):
                rows.append(irn)
                cols.append(ics[ik])
                data.append(nn[:,ik])

            ones = nm.ones( (nn.shape[0],), dtype = nm.float64 )
            for ik, il in enumerate(ir):
                rows.append(dim * ic + il)
                cols.append(ics[ik])
                data.append(ones)

        ## print rows
        ## print cols
        ## print data

        rows = nm.concatenate(rows)
        cols = nm.concatenate(cols)
        data = nm.concatenate(data)

        n_np_dof = n_nod * (dim - 1)
        mtx = sp.coo_matrix((data, (rows, cols)), shape=(n_nod * dim, n_np_dof))

        self.n_dof = n_np_dof
        self.mtx = mtx.tocsr()

        ## import pylab
        ## from sfepy.base.plotutils import spy
        ## spy( mtx )
        ## print mtx
        ## pylab.show()

class NormalDirectionOperator(LCBCOperator):
    """
    Transformation matrix operator for normal direction LCBCs.

    The substitution (in 3D) is:

    .. math::
        [u_1, u_2, u_3]^T = [n_1, n_2, n_3]^T w

    The new DOF is :math:`w`.
    """
    def __init__(self, name, nodes, region, field, dof_names, filename=None):
        Struct.__init__(self, name=name, nodes=nodes, dof_names=dof_names)

        dim = field.shape[0]
        assert_(len(dof_names) == dim)

        normals = compute_nodal_normals(nodes, region, field)

        if filename is not None:
            _save_normals(filename, normals, region, field.domain.mesh)

        n_nod, dim = normals.shape

        data = normals.ravel()
        rows = nm.arange(data.shape[0])
        cols = nm.repeat(nm.arange(n_nod), dim)

        mtx = sp.coo_matrix((data, (rows, cols)), shape=(n_nod * dim, n_nod))

        self.n_dof = n_nod
        self.mtx = mtx.tocsr()

class LCBCOperators(Container):
    """
    Container holding instances of LCBCOperator subclasses for a single
    variable.

    Parameters
    ----------
    name : str
        The object name.
    eq_map : EquationMap instance
        The equation mapping of the variable.
    offset : int
        The offset added to markers distinguishing the individual LCBCs.
    """
    def __init__(self, name, eq_map, offset):
        Container.__init__(self, name=name, eq_map=eq_map, offset=offset)

        self.eq_lcbc = nm.zeros((self.eq_map.n_eq,), dtype=nm.int32)
        self.markers = []
        self.n_transformed_dof = []
        self.n_op = 0
        self.ics = None

    def add_from_bc(self, bc, field):
        """
        Create a new LCBC operator described by `bc`, and add it to the
        container.

        Parameters
        ----------
        bc : LinearCombinationBC instance
            The LCBC condition description.
        field : Field instance
            The field of the variable.
        """
        region = bc.region
        dofs, kind = bc.dofs

        nmaster = region.get_field_nodes(field, merge=True)

        if kind == 'rigid':
            op = RigidOperator('%d_rigid' % len(self),
                               nmaster, field, dofs, self.eq_map.dof_names)

        elif kind == 'no_penetration':
            filename = get_default_attr(bc, 'filename', None)
            op = NoPenetrationOperator('%d_no_penetration' % len(self),
                                       nmaster, region, field, dofs,
                                       filename=filename)

        elif kind == 'normal_direction':
            filename = get_default_attr(bc, 'filename', None)
            op = NormalDirectionOperator('%d_normal_direction' % len(self),
                                         nmaster, region, field, dofs,
                                         filename=filename)

        self.append(op)

    def append(self, op):
        Container.append(self, op)

        eq = self.eq_map.eq
        meq = eq[expand_nodes_to_equations(op.nodes, op.dof_names,
                                           self.eq_map.dof_names)]
        assert_(nm.all(meq >= 0))

        op.treat_pbcs(meq)

        self.markers.append(self.offset + self.n_op + 1)
        self.eq_lcbc[meq] = self.markers[-1]

        self.n_transformed_dof.append(op.n_dof)
        self.n_op = len(self)

    def finalize(self):
        """
        Call this after all LCBCs of the variable have been added.

        Initializes the global column indices.
        """
        self.ics = nm.cumsum(nm.r_[0, self.n_transformed_dof])

def make_global_lcbc_operator(lcbc_ops, adi, new_only=False):
    """
    Assemble all LCBC operators into a single matrix.

    Returns
    -------
    mtx_lc : csr_matrix
        The global LCBC operator in the form of a CSR matrix.
    lcdi : DofInfo
        The global active LCBC-constrained DOF information.
    new_only : bool
        If True, the operator columns will contain only new DOFs.
    """
    n_dof = adi.ptr[-1]
    eq_lcbc = nm.zeros((n_dof,), dtype=nm.int32)

    n_dof_new = 0
    n_free = {}
    n_new = {}
    for var_name, lcbc_op in lcbc_ops.iteritems():
        ## print var_name, lcbc_op
        if lcbc_op is None: continue

        indx = adi.indx[var_name]
        eq_lcbc[indx] = lcbc_op.eq_lcbc

        n_free[var_name] = len(nm.where(lcbc_op.eq_lcbc == 0)[0])
        n_new[var_name] = nm.sum(lcbc_op.n_transformed_dof)

        n_dof_new += n_new[var_name]

    if n_dof_new == 0:
        return None, None

    ii = nm.nonzero( eq_lcbc )[0]
    n_constrained = ii.shape[0]
    n_dof_free = n_dof - n_constrained
    n_dof_reduced = n_dof_free + n_dof_new
    output( 'dofs: total %d, free %d, constrained %d, new %d'\
            % (n_dof, n_dof_free, n_constrained, n_dof_new) )
    output( ' -> reduced %d' % (n_dof_reduced) )

    lcdi = DofInfo('lcbc_active_state_dof_info')
    fdi = DofInfo('free_dof_info')
    ndi = DofInfo('new_dof_info')
    for var_name in adi.var_names:
        nf = n_free.get(var_name, adi.n_dof[var_name])
        nn = n_new.get(var_name, 0)
        fdi.append_raw(var_name, nf)
        ndi.append_raw(var_name, nn)
        lcdi.append_raw(var_name, nn + nf)

    assert_(lcdi.ptr[-1] == n_dof_reduced)

    rows = []
    cols = []
    data = []
    for var_name, lcbc_op in lcbc_ops.iteritems():
        if lcbc_op is None: continue

        if new_only:
            offset = ndi.indx[var_name].start

        else:
            offset = lcdi.indx[var_name].start + fdi.n_dof[var_name]

        for ii, op in enumerate(lcbc_op):
            indx = nm.where(eq_lcbc == lcbc_op.markers[ii])[0]
            icols = nm.arange(offset + lcbc_op.ics[ii],
                              offset + lcbc_op.ics[ii+1])

            if isinstance(op.mtx, sp.spmatrix):
                lr, lc, lv = sp.find(op.mtx)
                rows.append(indx[lr])
                cols.append(icols[lc])
                data.append(lv)

            else:
                irs, ics = nm.meshgrid(indx, icols)
                rows.append(irs.ravel())
                cols.append(ics.ravel())
                data.append(op.mtx.T.ravel())

    rows = nm.concatenate(rows)
    cols = nm.concatenate(cols)
    data = nm.concatenate(data)

    if new_only:
        mtx_lc = sp.coo_matrix((data, (rows, cols)),
                               shape=(n_dof, n_dof_new))

    else:
        mtx_lc = sp.coo_matrix((data, (rows, cols)),
                               shape=(n_dof, n_dof_reduced))

        ir = nm.where( eq_lcbc == 0 )[0]

        ic = nm.empty((n_dof_free,), dtype=nm.int32)
        for var_name in adi.var_names:
            ii = nm.arange(fdi.n_dof[var_name], dtype=nm.int32)
            ic[fdi.indx[var_name]] = lcdi.indx[var_name].start + ii

        mtx_lc2 = sp.coo_matrix((nm.ones((ir.shape[0],)), (ir, ic)),
                                shape=(n_dof, n_dof_reduced), dtype=nm.float64)


        mtx_lc = mtx_lc + mtx_lc2

    mtx_lc = mtx_lc.tocsr()

    ## import pylab
    ## from sfepy.base.plotutils import spy
    ## spy( mtx_lc )
    ## print mtx_lc
    ## pylab.show()

    return mtx_lc, lcdi

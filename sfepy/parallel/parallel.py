"""
Functions for a high-level PETSc-based parallelization.
"""
import numpy as nm

import sys, petsc4py
petsc4py.init(sys.argv)

from petsc4py import PETSc
from mpi4py import MPI

from sfepy.base.base import assert_, output, ordered_iteritems, Struct
from sfepy.discrete.common.region import Region
from sfepy.discrete.fem.fe_surface import FESurface

def get_inter_facets(domain, cell_tasks):
    """
    For each couple of neighboring task subdomains get the common boundary
    (interface) facets.
    """
    cmesh = domain.cmesh

    # Facet-to-cell connectivity.
    cmesh.setup_connectivity(cmesh.tdim - 1, cmesh.tdim)
    cfc = cmesh.get_conn(cmesh.tdim - 1, cmesh.tdim)

    # Facet tasks by cells in cfc.
    ftasks = cell_tasks[cfc.indices]

    # Mesh inner and surface facets.
    if_surf = cmesh.get_surface_facets()
    if_inner = nm.setdiff1d(nm.arange(cfc.num, dtype=nm.uint32), if_surf)

    # Facets in two tasks = inter-task region facets.
    if_inter = if_inner[nm.where(ftasks[cfc.offsets[if_inner]]
                                 != ftasks[cfc.offsets[if_inner] + 1])]
    aux = nm.c_[cfc.offsets[if_inter], cfc.offsets[if_inter] + 1]
    inter_tasks = ftasks[aux]

    # Fast version:
    # from sfepy.linalg import argsort_rows
    # ii = argsort_rows(inter_tasks)
    # inter_tasks[ii]
    # if_inter[ii]
    inter_facets = {}
    for ii, (i0, i1) in enumerate(inter_tasks):
        facet = if_inter[ii]

        ntasks = inter_facets.setdefault(i0, {})
        facets = ntasks.setdefault(i1, [])
        facets.append(facet)

        ntasks = inter_facets.setdefault(i1, {})
        facets = ntasks.setdefault(i0, [])
        facets.append(facet)

    return inter_facets

def create_task_dof_maps(field, cell_parts, inter_facets, is_overlap=True):
    """
    For each task list its inner and interface DOFs of the given field and
    create PETSc numbering that is consecutive in each subdomain.

    For each task, the DOF map has the following structure::

      [inner,
       [own_inter1, own_inter2, ...],
       [overlap_cells1, overlap_cells2, ...],
       n_task_total, task_offset]

    The overlapping cells are defined so that the system matrix corresponding
    to each task can be assembled independently, see [1]. TODO: Some "corner"
    cells may be added even if not needed - filter them out by using the PETSc
    DOFs range.

    [1] J. Sistek and F. Cirak. Parallel iterative solution of the
    incompressible Navier-Stokes equations with application to rotating
    wings. Submitted for publication, 2015
    """
    domain = field.domain
    cmesh = domain.cmesh

    id_map = nm.zeros(field.n_nod, dtype=nm.uint32)

    dof_maps = {}
    count = 0
    inter_count = 0
    ocs = nm.zeros(0, dtype=nm.int32)
    for ir, ntasks in ordered_iteritems(inter_facets):
        cregion = Region.from_cells(cell_parts[ir], domain, name='task_%d' % ir)
        domain.regions.append(cregion)
        dofs = field.get_dofs_in_region(cregion)

        rdof_map = dof_maps.setdefault(ir, [None, [], [], 0, 0])
        inter_dofs = []
        for ic, facets in ordered_iteritems(ntasks):
            cdof_map = dof_maps.setdefault(ic, [None, [], [], 0, 0])

            name = 'inter_%d_%d' % (ir, ic)
            ii = ir

            region = Region.from_facets(facets, domain, name,
                                        parent=cregion.name)
            region.update_shape()

            inter_dofs.append(field.get_dofs_in_region(region))

            ap = field.ap
            sd = FESurface('surface_data_%s' % region.name, region,
                           ap.efaces, ap.econn, field.region)
            econn = sd.get_connectivity()
            n_facet = econn.shape[0]

            ii2 = max(int(n_facet / 2), 1)

            dr = nm.unique(econn[:ii2])
            ii = nm.where((id_map[dr] == 0))[0]
            n_new = len(ii)
            if n_new:
                rdof_map[1].append(dr[ii])
                rdof_map[3] += n_new
                id_map[dr[ii]] = 1
                inter_count += n_new
                count += n_new

                if is_overlap:
                    ovs = cmesh.get_incident(0, region.facets[:ii2],
                                             cmesh.tdim - 1)
                    ocs = cmesh.get_incident(cmesh.tdim, ovs, 0)
                    rdof_map[2].append(ocs.astype(nm.int32))

            dc = nm.unique(econn[ii2:])
            ii = nm.where((id_map[dc] == 0))[0]
            n_new = len(ii)
            if n_new:
                cdof_map[1].append(dc[ii])
                cdof_map[3] += n_new
                id_map[dc[ii]] = 1
                inter_count += n_new
                count += n_new

                if is_overlap:
                    ovs = cmesh.get_incident(0, region.facets[ii2:],
                                             cmesh.tdim - 1)
                    ocs = cmesh.get_incident(cmesh.tdim, ovs, 0)
                    cdof_map[2].append(ocs.astype(nm.int32))

        domain.regions.pop() # Remove the cell region.

        inner_dofs = nm.setdiff1d(dofs, nm.concatenate(inter_dofs))
        n_inner = len(inner_dofs)
        rdof_map[3] += n_inner
        assert_(nm.all(id_map[inner_dofs] == 0))
        id_map[inner_dofs] = 1
        count += n_inner

        rdof_map[0] = inner_dofs

    offset = 0
    overlap_cells = []
    for ir, dof_map in ordered_iteritems(dof_maps):
        n_owned = dof_map[3]
        output(n_owned, offset)

        i0 = len(dof_map[0])
        id_map[dof_map[0]] = nm.arange(offset, offset + i0, dtype=nm.uint32)
        for aux in dof_map[1]:
            i1 = len(aux)
            id_map[aux] = nm.arange(offset + i0, offset + i0 + i1,
                                    dtype=nm.uint32)
            i0 += i1

        if len(dof_map[2]):
            ocs = nm.unique(nm.concatenate(dof_map[2]))

        else:
            ocs = nm.zeros(0, dtype=nm.int32)

        overlap_cells.append(ocs)

        assert_(i0 == n_owned)

        dof_map[4] = offset
        offset += n_owned

    if not len(overlap_cells):
        overlap_cells.append( nm.zeros(0, dtype=nm.int32))

    return dof_maps, id_map, overlap_cells

def distribute_field_dofs(field, cell_parts, cell_tasks, is_overlap=True,
                          comm=None, verbose=False):
    """
    Distribute the owned cells and DOFs of the given field to all tasks.

    The DOFs use the PETSc ordering and are in form of a connectivity, so that
    each task can easily identify them with the DOFs of the original global
    ordering or local ordering.
    """
    if comm is None:
        comm = PETSc.COMM_WORLD

    size = comm.size
    mpi = comm.tompi4py()

    if comm.rank == 0:
        inter_facets = get_inter_facets(field.domain, cell_tasks)

        aux = create_task_dof_maps(field, cell_parts, inter_facets,
                                   is_overlap=is_overlap)
        dof_maps, id_map, overlap_cells = aux

        n_cell_parts = [len(ii) for ii in cell_parts]
        output('numbers of cells in tasks (without overlaps):',
               n_cell_parts, verbose=verbose)
        assert_(sum(n_cell_parts) == field.domain.mesh.n_el)
        assert_(nm.all(n_cell_parts > 0))

        # Send subdomain data to other tasks.
        for it in xrange(1, size):
            # Send owned and overlap cells.
            cells = nm.union1d(cell_parts[it], overlap_cells[it])
            mpi.send(len(cells), it)
            mpi.Send([cells, MPI.INTEGER4], it)

            dof_map = dof_maps[it]

            # Send owned petsc_dofs range.
            mpi.send(dof_map[4], it)
            mpi.send(dof_map[4] + dof_map[3], it)

            # Send petsc_dofs of global_dofs.
            global_dofs = field.ap.econn[cells]
            petsc_dofs_conn = id_map[global_dofs]
            mpi.send(petsc_dofs_conn.shape[0], it)
            mpi.send(petsc_dofs_conn.shape[1], it)
            mpi.Send([petsc_dofs_conn, MPI.INTEGER4], it)

        cells = nm.union1d(cell_parts[0], overlap_cells[0])
        n_cell = len(cells)

        global_dofs = field.ap.econn[cells]

        if 0 in dof_maps:
            dof_map = dof_maps[0]
            petsc_dofs_range = (dof_map[4], dof_map[4] + dof_map[3])
            petsc_dofs_conn = id_map[global_dofs]

        else:
            petsc_dofs_range = (0, global_dofs.max() + 1)
            petsc_dofs_conn = global_dofs

    else:
        # Receive owned cells.
        n_cell = mpi.recv(source=0)
        cells = nm.empty(n_cell, dtype=nm.int32)
        mpi.Recv([cells, MPI.INTEGER4], source=0)

        # Receive owned petsc_dofs range.
        i0 = mpi.recv(source=0)
        i1 = mpi.recv(source=0)
        petsc_dofs_range = (i0, i1)

        # Receive petsc_dofs of global_dofs.
        n_cell = mpi.recv(source=0)
        n_cdof = mpi.recv(source=0)
        petsc_dofs_conn = nm.empty((n_cell, n_cdof), dtype=nm.int32)
        mpi.Recv([petsc_dofs_conn, MPI.INTEGER4], source=0)

        dof_maps = id_map = None

    if verbose:
        output('n_cell:', n_cell)
        output('cells:', cells)
        output('owned petsc DOF range:', petsc_dofs_range,
               petsc_dofs_range[1] - petsc_dofs_range[0])
        aux = nm.unique(petsc_dofs_conn)
        output('local petsc DOFs (owned + shared):', aux, len(aux))

    return cells, petsc_dofs_range, petsc_dofs_conn, dof_maps, id_map

def get_local_ordering(field_i, petsc_dofs_conn):
    """
    Get PETSc DOFs in the order of local DOFs of the localized field `field_i`.
    """
    petsc_dofs = nm.empty(field_i.n_nod, dtype=nm.int32)
    econn = field_i.ap.econn
    petsc_dofs[econn] = petsc_dofs_conn

    return petsc_dofs

def get_sizes(petsc_dofs_range, n_dof, n_components):
    """
    Get (local, total) sizes of a vector and local equation range.
    """
    drange = tuple(n_components * nm.asarray(petsc_dofs_range))
    n_loc = drange[1] - drange[0]
    n_all_dof = n_dof * n_components
    sizes = (n_loc, n_all_dof)

    return sizes, drange

def expand_dofs(dofs, n_components):
    """
    Expand DOFs to equation numbers.
    """
    edofs = nm.empty(n_components * dofs.shape[0], nm.int32)
    for idof in xrange(n_components):
        aux = n_components * dofs + idof
        edofs[idof::n_components] = aux

    return edofs

def create_prealloc_data(mtx, pdofs, drange, verbose=False):
    """
    Create CSR preallocation data for a PETSc matrix based on the owned PETSc
    DOFs and a local matrix with EBCs not applied.
    """
    owned_dofs = nm.where((pdofs >= drange[0]) & (pdofs < drange[1]))[0]
    owned_dofs = owned_dofs.astype(nm.int32)
    output('owned local DOFs:', owned_dofs, verbose=verbose)

    ii = nm.argsort(pdofs[owned_dofs])
    aux = mtx[owned_dofs[ii]]
    mtx_prealloc = Struct(indptr=aux.indptr, indices=pdofs[aux.indices])

    return mtx_prealloc

def create_petsc_matrix(sizes, mtx_prealloc=None, comm=None):
    """
    Create and allocate a PETSc matrix.
    """
    if comm is None:
        comm = PETSc.COMM_WORLD

    pmtx = PETSc.Mat()
    pmtx.create(comm)
    pmtx.setType('aij')

    pmtx.setSizes((sizes, sizes))

    if mtx_prealloc is not None:
        pmtx.setPreallocationCSR((mtx_prealloc.indptr, mtx_prealloc.indices))

    pmtx.setUp()

    return pmtx

def apply_ebc_to_matrix(mtx, eq_map):
    """
    Apply to matrix rows: zeros to non-diagonal entries, one to the diagonal.
    """
    ebc_rows = eq_map.eq_ebc

    data, prows, cols = mtx.data, mtx.indptr, mtx.indices
    for ir in ebc_rows:
        for ic in xrange(prows[ir], prows[ir + 1]):
            if (cols[ic] == ir):
                data[ic] = 1.0

            else:
                data[ic] = 0.0

def assemble_to_petsc(pmtx, prhs, mtx, rhs, pdofs, drange, is_overlap=True,
                      comm=None, verbose=False):
    """
    Assemble local CSR matrix and right-hand side vector to PETSc counterparts.

    WIP
    ---
    Try Mat.setValuesCSR() - no lgmap - filtering vectorized?
    """
    import time

    if comm is None:
        comm = PETSc.COMM_WORLD

    lgmap = PETSc.LGMap().create(pdofs, comm=comm)

    if is_overlap:
        pmtx.setLGMap(lgmap, lgmap)

        data, prows, cols = mtx.data, mtx.indptr, mtx.indices

        output('setting matrix values...', verbose=verbose)
        tt = time.clock()
        for ir, rdof in enumerate(pdofs):
            if (rdof < drange[0]) or (rdof >= drange[1]): continue

            for ic in xrange(prows[ir], prows[ir + 1]):
                # output(ir, rdof, cols[ic])
                pmtx.setValueLocal(ir, cols[ic], data[ic],
                                   PETSc.InsertMode.INSERT_VALUES)
        output('...done in', time.clock() - tt, verbose=verbose)

        output('assembling matrix...', verbose=verbose)
        tt = time.clock()
        pmtx.assemble()
        output('...done in', time.clock() - tt, verbose=verbose)

        prhs.setLGMap(lgmap)
        output('setting rhs values...', verbose=verbose)
        tt = time.clock()
        for ir, rdof in enumerate(pdofs):
            if (rdof < drange[0]) or (rdof >= drange[1]): continue
            prhs.setValueLocal(ir, rhs[ir],
                               PETSc.InsertMode.INSERT_VALUES)
        output('...done in', time.clock() - tt, verbose=verbose)

        output('assembling rhs...', verbose=verbose)
        tt = time.clock()
        prhs.assemble()
        output('...done in', time.clock() - tt, verbose=verbose)

    else:
        pmtx.setLGMap(lgmap, lgmap)
        output('setting matrix values...', verbose=verbose)
        tt = time.clock()
        pmtx.setValuesLocalCSR(mtx.indptr, mtx.indices, mtx.data,
                               PETSc.InsertMode.ADD_VALUES)
        output('...done in', time.clock() - tt, verbose=verbose)

        output('assembling matrix...', verbose=verbose)
        tt = time.clock()
        pmtx.assemble()
        output('...done in', time.clock() - tt, verbose=verbose)

        prhs.setLGMap(lgmap)
        output('setting rhs values...', verbose=verbose)
        tt = time.clock()
        prhs.setValuesLocal(nm.arange(len(rhs), dtype=nm.int32), rhs,
                            PETSc.InsertMode.ADD_VALUES)
        output('...done in', time.clock() - tt, verbose=verbose)

        output('assembling rhs...', verbose=verbose)
        tt = time.clock()
        prhs.assemble()
        output('...done in', time.clock() - tt, verbose=verbose)

"""
Functions for a high-level PETSc-based parallelization.
"""
from __future__ import absolute_import
import os

import numpy as nm
from six.moves import range

def init_petsc_args():
    try:
        import sys, petsc4py

    except ImportError:
        return

    argv = [arg for arg in sys.argv if arg not in ['-h', '--help']]
    petsc4py.init(argv)

init_petsc_args()

try:
    from petsc4py import PETSc

except (ModuleNotFoundError, ImportError):
    PETSc = None

try:
    from mpi4py import MPI

except (ModuleNotFoundError, ImportError):
    MPI = None

from sfepy.base.base import assert_, output, ordered_iteritems, Struct
from sfepy.base.timing import Timer
from sfepy.discrete.common.region import Region
from sfepy.discrete.fem.fe_surface import FESurface

def partition_mesh(mesh, n_parts, use_metis=True, verbose=False):
    """
    Partition the mesh cells into `n_parts` subdomains, using metis, if
    available.
    """
    output('partitioning mesh into %d subdomains...' % n_parts, verbose=verbose)
    timer = Timer(start=True)

    if use_metis:
        try:
            from pymetis import part_graph

        except ImportError:
            output('pymetis is not available, using naive partitioning!')
            part_graph = None

    if use_metis and (part_graph is not None):
        cmesh = mesh.cmesh
        cmesh.setup_connectivity(cmesh.dim, cmesh.dim)
        graph = cmesh.get_conn(cmesh.dim, cmesh.dim)

        cuts, cell_tasks = part_graph(n_parts, xadj=graph.offsets.astype(int),
                                      adjncy=graph.indices.astype(int))
        cell_tasks = nm.array(cell_tasks, dtype=nm.int32)

    else:
        ii = nm.arange(n_parts)
        n_cell_parts = mesh.n_el // n_parts + ((mesh.n_el % n_parts) > ii)
        output('cell counts:', n_cell_parts, verbose=verbose)
        assert_(sum(n_cell_parts) == mesh.n_el)
        assert_(nm.all(n_cell_parts > 0))

        offs = nm.cumsum(nm.r_[0, n_cell_parts])
        cell_tasks = nm.digitize(nm.arange(offs[-1]), offs) - 1

    output('...done in', timer.stop(), verbose=verbose)

    return cell_tasks

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

def create_task_dof_maps(field, cell_tasks, inter_facets, is_overlap=True,
                         use_expand_dofs=False, save_inter_regions=False,
                         output_dir=None):
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

    When debugging domain partitioning problems, it is advisable to set
    `save_inter_regions` to True to save the task interfaces as meshes as well
    as vertex-based markers - to be used only with moderate problems and small
    numbers of tasks.

    [1] J. Sistek and F. Cirak. Parallel iterative solution of the
    incompressible Navier-Stokes equations with application to rotating
    wings. Submitted for publication, 2015
    """
    domain = field.domain
    cmesh = domain.cmesh

    if use_expand_dofs:
        id_map = nm.zeros(field.n_nod * field.n_components, dtype=nm.uint32)

    else:
        id_map = nm.zeros(field.n_nod, dtype=nm.uint32)

    def _get_dofs_region(field, region):
        dofs = field.get_dofs_in_region(region)
        if use_expand_dofs:
            dofs = expand_dofs(dofs, field.n_components)
        return dofs

    def _get_dofs_conn(field, conn):
        dofs = nm.unique(conn)
        if use_expand_dofs:
            dofs = expand_dofs(dofs, field.n_components)
        return dofs

    dof_maps = {}
    count = 0
    inter_count = 0
    ocs = nm.zeros(0, dtype=nm.int32)
    cell_parts = []
    for ir, ntasks in ordered_iteritems(inter_facets):
        cells = nm.where(cell_tasks == ir)[0].astype(nm.int32)
        cell_parts.append(cells)

        cregion = Region.from_cells(cells, domain, name='task_%d' % ir)
        domain.regions.append(cregion)
        dofs = _get_dofs_region(field, cregion)

        rdof_map = dof_maps.setdefault(ir, [None, [], [], 0, 0])
        inter_dofs = []
        for ic, facets in ordered_iteritems(ntasks):
            cdof_map = dof_maps.setdefault(ic, [None, [], [], 0, 0])

            name = 'inter_%d_%d' % (ir, ic)
            ii = ir

            region = Region.from_facets(facets, domain, name,
                                        parent=cregion.name)
            region.update_shape()

            if save_inter_regions:
                output_dir = output_dir if output_dir is not None else '.'
                filename = os.path.join(output_dir, '%s.mesh' % name)

                aux = domain.mesh.from_region(region, domain.mesh,
                                              is_surface=True)
                aux.write(filename, io='auto')

                mask = nm.zeros((domain.mesh.n_nod, 1), dtype=nm.float64)
                mask[region.vertices] = 1
                out = {name : Struct(name=name, mode='vertex', data=mask)}
                filename = os.path.join(output_dir, '%s.h5' % name)
                domain.mesh.write(filename, out=out, io='auto')

            inter_dofs.append(_get_dofs_region(field, region))

            sd = FESurface('surface_data_%s' % region.name, region,
                           field.efaces, field.econn, field.region)
            econn = sd.get_connectivity()
            n_facet = econn.shape[0]

            ii2 = max(int(n_facet / 2), 1)

            dr = _get_dofs_conn(field, econn[:ii2])
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

            dc = _get_dofs_conn(field, econn[ii2:])
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

    if not len(dof_maps):
        dofs = _get_dofs_region(field, field.region)
        dof_maps[0] = [dofs, [], [], len(dofs), 0]
        id_map[:] = nm.arange(len(dofs), dtype=nm.uint32)

    if not len(cell_parts):
        cell_parts.append(nm.arange(len(cell_tasks), dtype=nm.int32))

    if not len(overlap_cells):
        overlap_cells.append(nm.zeros(0, dtype=nm.int32))

    return dof_maps, id_map, cell_parts, overlap_cells

def verify_task_dof_maps(dof_maps, id_map, field, use_expand_dofs=False,
                         verbose=False):
    """
    Verify the counts and values of DOFs in `dof_maps` and `id_map`
    corresponding to `field`.

    Returns the vector with a task number for each DOF.
    """
    timer = Timer(start=True)
    if verbose:
        output('verifying...')
        output('total number of DOFs:', field.n_nod)
        output('number of tasks:', len(dof_maps))

    count = count2 = 0
    dofs = []
    if use_expand_dofs:
        vec = nm.empty(field.n_nod * field.n_components, dtype=nm.float64)

    else:
        vec = nm.empty(field.n_nod, dtype=nm.float64)
    for ir, dof_map in ordered_iteritems(dof_maps):
        n_owned = dof_map[3]
        offset = dof_map[4]
        o2 = offset + n_owned

        if verbose:
            output('task %d: %d owned on offset %d' % (ir, n_owned, offset))

        if not use_expand_dofs:
            aux = dof_map[0]
            assert_(nm.all((id_map[aux] >= offset) & (id_map[aux] < o2)))

        count2 += dof_map[3]

        count += len(dof_map[0])

        dofs.append(dof_map[0])
        vec[dof_map[0]] = ir
        for aux in dof_map[1]:
            if not use_expand_dofs:
                assert_(nm.all((id_map[aux] >= offset) & (id_map[aux] < o2)))

            count += len(aux)
            dofs.append(aux)
            vec[aux] = ir

    dofs = nm.concatenate(dofs)

    n_dof = vec.shape[0]

    assert_(n_dof == len(dofs))
    if not expand_dofs:
        assert_(nm.all(nm.sort(dofs) == nm.sort(id_map)))

    dofs = nm.unique(dofs)
    assert_(n_dof == len(dofs))

    assert_(n_dof == dofs[-1] + 1)
    assert_(n_dof == count)
    assert_(n_dof == count2)
    assert_(n_dof == len(id_map))
    assert_(n_dof == len(nm.unique(id_map)))

    output('...done in', timer.stop(), verbose=verbose)

    return vec

def distribute_field_dofs(field, gfd, use_expand_dofs=False,
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
        dof_maps = gfd.dof_maps
        id_map = gfd.id_map

        # Send subdomain data to other tasks.
        for it in range(1, size):
            # Send owned and overlap cells.
            cells = nm.union1d(gfd.cell_parts[it], gfd.overlap_cells[it])
            mpi.send(len(cells), it)
            mpi.Send([cells, MPI.INTEGER4], it)

            dof_map = dof_maps[it]

            # Send owned petsc_dofs range.
            mpi.send(gfd.coffsets[it], it)
            mpi.send(gfd.coffsets[it] + dof_map[3], it)

            # Send petsc_dofs of global_dofs.
            global_dofs = field.econn[cells]
            if use_expand_dofs:
                global_dofs = expand_dofs(global_dofs, field.n_components)
            petsc_dofs_conn = id_map[global_dofs]
            mpi.send(petsc_dofs_conn.shape[0], it)
            mpi.send(petsc_dofs_conn.shape[1], it)
            mpi.Send([petsc_dofs_conn, MPI.INTEGER4], it)

        cells = nm.union1d(gfd.cell_parts[0], gfd.overlap_cells[0])
        n_cell = len(cells)

        global_dofs = field.econn[cells]
        if use_expand_dofs:
            global_dofs = expand_dofs(global_dofs, field.n_components)

        dof_map = dof_maps[0]
        petsc_dofs_range = (gfd.coffsets[0], gfd.coffsets[0] + dof_map[3])
        petsc_dofs_conn = id_map[global_dofs]

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
        output('field %s:' % field.name)
        output('n_cell:', n_cell)
        output('cells:', cells)
        output('owned petsc DOF range:', petsc_dofs_range,
               petsc_dofs_range[1] - petsc_dofs_range[0])
        aux = nm.unique(petsc_dofs_conn)
        output('%d local petsc DOFs (owned + shared):' % len(aux), aux)

    return cells, petsc_dofs_range, petsc_dofs_conn, dof_maps, id_map

def distribute_fields_dofs(fields, cell_tasks, is_overlap=True,
                           use_expand_dofs=False, save_inter_regions=False,
                           output_dir=None, comm=None, verbose=False):
    """
    Distribute the owned cells and DOFs of the given field to all tasks.

    Uses interleaved PETSc numbering in each task, i.e., the PETSc DOFs of each
    tasks are consecutive and correspond to the first field DOFs block followed
    by the second etc.

    Expand DOFs to equations if `use_expand_dofs` is True.
    """
    if comm is None:
        comm = PETSc.COMM_WORLD

    size = comm.size

    if comm.rank == 0:
        gfds = []
        inter_facets = get_inter_facets(fields[0].domain, cell_tasks)
        for field in fields:
            aux = create_task_dof_maps(field, cell_tasks, inter_facets,
                                       is_overlap=is_overlap,
                                       use_expand_dofs=use_expand_dofs,
                                       save_inter_regions=save_inter_regions,
                                       output_dir=output_dir)
            cell_parts = aux[2]
            n_cell_parts = [len(ii) for ii in cell_parts]
            output('numbers of cells in tasks (without overlaps):',
                   n_cell_parts, verbose=verbose)
            assert_(sum(n_cell_parts) == field.domain.mesh.n_el)
            assert_(nm.all(nm.array(n_cell_parts) > 0))

            gfd = Struct(name='global field %s distribution' % field.name,
                         dof_maps=aux[0], id_map=aux[1],
                         cell_parts=aux[2], overlap_cells=aux[3],
                         coffsets=nm.empty(size, dtype=nm.int32))
            gfds.append(gfd)

        # Initialize composite offsets of DOFs.
        if len(fields) > 1:
            # Renumber id_maps for field inter-leaving.
            offset = 0
            for ir in range(size):
                for ii, gfd in enumerate(gfds):
                    dof_map = gfd.dof_maps[ir]
                    n_owned = dof_map[3]
                    off = dof_map[4]

                    iown = nm.concatenate([dof_map[0]] + dof_map[1])
                    gfd.id_map[iown] += offset - off
                    gfd.coffsets[ir] = offset

                    offset += n_owned

        else:
            gfd = gfds[0]
            gfd.coffsets[:] = [gfd.dof_maps[ir][4] for ir in range(size)]

    else:
        gfds = [None] * len(fields)

    lfds = []
    for ii, field in enumerate(fields):
        aux = distribute_field_dofs(field, gfds[ii],
                                    use_expand_dofs=use_expand_dofs,
                                    comm=comm, verbose=verbose)
        lfd = Struct(name='local field %s distribution' % field.name,
                      cells=aux[0], petsc_dofs_range=aux[1],
                      petsc_dofs_conn=aux[2])
        lfds.append(lfd)

    return lfds, gfds

def get_local_ordering(field_i, petsc_dofs_conn, use_expand_dofs=False):
    """
    Get PETSc DOFs in the order of local DOFs of the localized field `field_i`.

    Expand DOFs to equations if `use_expand_dofs` is True.
    """
    if use_expand_dofs:
        petsc_dofs = nm.empty(field_i.n_nod * field_i.n_components,
                              dtype=nm.int32)

    else:
        petsc_dofs = nm.empty(field_i.n_nod, dtype=nm.int32)
    econn = field_i.econn
    if use_expand_dofs:
        econn = expand_dofs(econn, field_i.n_components)
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

def get_composite_sizes(lfds):
    """
    Get (local, total) sizes of a vector and local equation range for a
    composite matrix built from field blocks described by `lfds` local field
    distributions information.
    """
    sizes = tuple(sum(ii) for ii in zip(*[ii.sizes for ii in lfds]))
    drange = (lfds[0].drange[0], lfds[-1].drange[1])

    return sizes, drange

def setup_composite_dofs(lfds, fields, local_variables, verbose=False):
    """
    Setup composite DOFs built from field blocks described by `lfds` local
    field distributions information.

    Returns (local, total) sizes of a vector, local equation range for a
    composite matrix, and the local ordering of composite PETSc DOFs,
    corresponding to `local_variables` (must be in the order of `fields`!).
    """
    for ii, variable in enumerate(local_variables.iter_state(ordered=True)):
        output('field %s:' % fields[ii].name, verbose=verbose)
        lfd = lfds[ii]
        output('PETSc DOFs range:', lfd.petsc_dofs_range, verbose=verbose)

        n_cdof = fields[ii].n_nod * fields[ii].n_components
        lfd.sizes, lfd.drange = get_sizes(lfd.petsc_dofs_range, n_cdof, 1)
        output('sizes, drange:', lfd.sizes, lfd.drange, verbose=verbose)

        lfd.petsc_dofs = get_local_ordering(variable.field,
                                            lfd.petsc_dofs_conn,
                                            use_expand_dofs=True)
        output('%d petsc dofs:' % len(lfd.petsc_dofs), lfd.petsc_dofs,
               verbose=verbose)

    sizes, drange = get_composite_sizes(lfds)
    output('composite sizes:', sizes, verbose=verbose)
    output('composite drange:', drange, verbose=verbose)

    pdofs = nm.concatenate([ii.petsc_dofs for ii in lfds])
    output('%d composite pdofs:' % len(pdofs), pdofs, verbose=verbose)

    return sizes, drange, pdofs

def expand_dofs(dofs, n_components):
    """
    Expand DOFs to equation numbers.
    """
    if dofs.ndim > 1:
        sh = dofs.shape
        dofs = dofs.ravel()

    else:
        sh = None

    edofs = nm.empty(n_components * dofs.shape[0], nm.int32)
    for idof in range(n_components):
        aux = n_components * dofs + idof
        edofs[idof::n_components] = aux

    if sh is not None:
        edofs.shape = sh[:-1] + (-1,)

    return edofs

def call_in_rank_order(fun, comm=None):
    """
    Call a function `fun` task by task in the task rank order.
    """
    if comm is None:
        comm = PETSc.COMM_WORLD

    for rank in range(comm.size):
        if rank == comm.rank:
            fun(rank, comm)
        comm.barrier()

def view_petsc_local(data, name='data', viewer=None, comm=None):
    """
    View local PETSc `data` called `name`. The data object has to have
    `.view()` method.
    """
    def _view(rank, comm):
        output('contents of local %s on rank %d:' % (name, rank))
        data.view(viewer=viewer)

    call_in_rank_order(_view, comm=comm)

def create_local_petsc_vector(pdofs):
    """
    Create a local PETSc vector with the size corresponding to `pdofs`.
    """
    pvec_i = PETSc.Vec().create(comm=PETSc.COMM_SELF)
    pvec_i.setSizes(len(pdofs))
    pvec_i.setUp()

    return pvec_i

def create_gather_scatter(pdofs, pvec_i, pvec, comm=None):
    """
    Create the ``gather()`` function for updating a global PETSc vector from
    local ones and the ``scatter()`` function for updating local PETSc vectors
    from the global one.
    """
    if comm is None:
        comm = PETSc.COMM_WORLD

    isg = PETSc.IS().createGeneral(pdofs, comm=comm)
    g2l = PETSc.Scatter().create(pvec, isg, pvec_i, None)

    def scatter(pvec_i, pvec):
        """
        Scatter `pvec` to `pvec_i`.
        """
        g2l.scatter(pvec, pvec_i, PETSc.InsertMode.INSERT)

    def gather(pvec, pvec_i):
        """
        Gather `pvec_i` to `pvec`.
        """
        g2l.scatter(pvec_i, pvec, PETSc.InsertMode.INSERT,
                    PETSc.ScatterMode.REVERSE)

    return gather, scatter

def create_gather_to_zero(pvec):
    """
    Create the ``gather_to_zero()`` function for collecting the global PETSc
    vector on the task of rank zero.
    """
    g20, pvec_full =  PETSc.Scatter().toZero(pvec)

    def gather_to_zero(pvec):
        """
        Return the global PETSc vector, corresponding to `pvec`, on the task of
        rank zero. The vector is reused between calls!
        """
        g20.scatter(pvec, pvec_full, PETSc.InsertMode.INSERT,
                    PETSc.ScatterMode.FORWARD)

        return pvec_full

    return gather_to_zero

def create_prealloc_data(mtx, pdofs, drange, verbose=False):
    """
    Create CSR preallocation data for a PETSc matrix based on the owned PETSc
    DOFs and a local matrix with EBCs not applied.
    """
    owned_dofs = nm.where((pdofs >= drange[0]) & (pdofs < drange[1]))[0]
    owned_dofs = owned_dofs.astype(nm.int32)
    output('%d owned local DOFs:' % len(owned_dofs), owned_dofs,
           verbose=verbose)

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

def create_petsc_system(mtx, sizes, pdofs, drange, is_overlap=True,
                        comm=None, verbose=False):
    """
    Create and pre-allocate (if `is_overlap` is True) a PETSc matrix and
    related solution and right-hand side vectors.
    """
    if comm is None:
        comm = PETSc.COMM_WORLD

    if is_overlap:
        mtx.data[:] = 1
        mtx_prealloc = create_prealloc_data(mtx, pdofs, drange,
                                            verbose=True)
        pmtx = create_petsc_matrix(sizes, mtx_prealloc, comm=comm)

    else:
        pmtx = create_petsc_matrix(sizes, comm=comm)

    own_range = pmtx.getOwnershipRange()
    output('pmtx ownership:', own_range, verbose=verbose)
    assert_(own_range == drange)

    psol, prhs = pmtx.getVecs()

    own_range = prhs.getOwnershipRange()
    output('prhs ownership:', own_range, verbose=verbose)
    assert_(own_range == drange)

    return pmtx, psol, prhs

def assemble_rhs_to_petsc(prhs, rhs, pdofs, drange, is_overlap=True,
                          comm=None, verbose=False):
    """
    Assemble a local right-hand side vector to a global PETSc vector.
    """
    if comm is None:
        comm = PETSc.COMM_WORLD

    timer = Timer()
    if is_overlap:
        output('setting rhs values...', verbose=verbose)
        timer.start()
        rdofs = nm.where((pdofs < drange[0]) | (pdofs >= drange[1]), -1, pdofs)
        prhs.setOption(prhs.Option.IGNORE_NEGATIVE_INDICES, True)
        prhs.setValues(rdofs, rhs, PETSc.InsertMode.INSERT_VALUES)
        output('...done in', timer.stop(), verbose=verbose)

        output('assembling rhs...', verbose=verbose)
        timer.start()
        prhs.assemble()
        output('...done in', timer.stop(), verbose=verbose)

    else:
        output('setting rhs values...', verbose=verbose)
        timer.start()
        prhs.setValues(pdofs, rhs, PETSc.InsertMode.ADD_VALUES)
        output('...done in', timer.stop(), verbose=verbose)

        output('assembling rhs...', verbose=verbose)
        timer.start()
        prhs.assemble()
        output('...done in', timer.stop(), verbose=verbose)

def assemble_mtx_to_petsc(pmtx, mtx, pdofs, drange, is_overlap=True,
                          comm=None, verbose=False):
    """
    Assemble a local CSR matrix to a global PETSc matrix.
    """
    if comm is None:
        comm = PETSc.COMM_WORLD

    timer = Timer()

    lgmap = PETSc.LGMap().create(pdofs, comm=comm)
    pmtx.setLGMap(lgmap, lgmap)
    if is_overlap:
        output('setting matrix values...', verbose=verbose)
        timer.start()
        mask = (pdofs < drange[0]) | (pdofs >= drange[1])
        nnz_per_row = nm.diff(mtx.indptr)
        mtx2 = mtx.copy()
        mtx2.data[nm.repeat(mask, nnz_per_row)] = 0
        mtx2.eliminate_zeros()
        pmtx.setValuesLocalCSR(mtx2.indptr, mtx2.indices, mtx2.data,
                               PETSc.InsertMode.INSERT_VALUES)
        output('...done in', timer.stop(), verbose=verbose)

        output('assembling matrix...', verbose=verbose)
        timer.start()
        pmtx.assemble()
        output('...done in', timer.stop(), verbose=verbose)


    else:
        output('setting matrix values...', verbose=verbose)
        timer.start()
        pmtx.setValuesLocalCSR(mtx.indptr, mtx.indices, mtx.data,
                               PETSc.InsertMode.ADD_VALUES)
        output('...done in', timer.stop(), verbose=verbose)

        output('assembling matrix...', verbose=verbose)
        timer.start()
        pmtx.assemble()
        output('...done in', timer.stop(), verbose=verbose)

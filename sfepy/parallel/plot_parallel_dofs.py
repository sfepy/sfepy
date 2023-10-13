"""
Functions to visualize the partitioning of a domain and a field DOFs.
"""
import os

import numpy as nm

from sfepy.base.base import ordered_iteritems
import sfepy.postprocess.plot_cmesh as pc
import sfepy.postprocess.plot_dofs as pd

def mark_subdomains(ax, cmesh, cell_tasks,
                    size=None, icolor=0, alpha=1.0, mask=False):
    """
    Mark cells of subdomains corresponding to each task by a different color.
    Plots nothing in 3D.
    """
    if size is None:
        size = cell_tasks.max() + 1

    coors = cmesh.coors
    dim = cmesh.dim

    ax = pd._get_axes(ax, dim)

    if dim == 3:
        return ax

    conn = cmesh.get_conn(dim, 0)

    color = nm.zeros(4)
    color[-1] = alpha
    for ic, vertices in enumerate(conn.indices.reshape((conn.num, -1))):
        cv = coors[vertices]
        if dim == 2:
            if mask:
                cc = cv.mean(0)
                cv = cc + 0.3 * (cv - cc)

            if not mask or cell_tasks[ic] > 0:
                color[icolor] = (float(cell_tasks[ic]) + 1) / size
                ax.fill(*cv.T, color=color)

    return ax

def label_dofs(ax, coors, dofs, colors):
    """
    Label DOFs using the given colors.
    """
    from sfepy.postprocess.plot_dofs import _get_axes

    dim = coors.shape[1]
    ax = _get_axes(ax, dim)

    for gdof in dofs:
        cd = coors[gdof]
        ax.text(*cd.T, s='%d' % gdof,
                color=colors[gdof], fontsize=12, weight='bold')

    return ax

def plot_partitioning(axs, field, cell_tasks, gfd, output_dir, size):
    """
    Plot the partitioning of the domain and field DOFs.
    """
    mesh = field.domain.mesh

    ax = pc.plot_wireframe(axs[0], mesh.cmesh)
    coors = field.get_coor()
    econn = field.econn
    ax = pd.plot_global_dofs(ax, coors, econn)
    ax.set_title('global DOFs')
    ax.figure.savefig(os.path.join(output_dir, 'global_dofs.png'),
                      bbox_inches='tight')

    ax = pc.plot_wireframe(axs[1], mesh.cmesh)
    fig = ax.figure
    coors = field.get_coor()
    econn = field.econn

    id_map = gfd.id_map

    colors = nm.zeros((field.n_nod, 4))
    for ir, dof_map in ordered_iteritems(gfd.dof_maps):
        aux = id_map[dof_map[0]]
        colors[aux] = [0, 0, float(ir + 1) / size, 0.6]
        for aux in dof_map[1]:
            colors[id_map[aux]] = [0, 0, float(ir + 1) / size, 0.9]

    from sfepy.discrete.fem.utils import prepare_translate
    aux = prepare_translate(id_map[econn], econn)
    ax = label_dofs(ax, coors[aux], id_map, colors)

    mark_subdomains(ax, mesh.cmesh, cell_tasks, size, 0, 0.7)

    ax.set_title('subdomains of tasks and PETSc DOFs')

    fig.savefig(os.path.join(output_dir, 'petsc_dofs.png'),
                bbox_inches='tight')

    ax.set_title('')

    axis = ax.axis()
    for ir, ocells in enumerate(gfd.overlap_cells):
        aux = nm.zeros_like(cell_tasks)
        aux[ocells] = 10
        aux[gfd.cell_parts[ir]] = 1

        ax = fig.add_axes(ax.get_position(), frameon=False, label='aux')
        mark_subdomains(ax, mesh.cmesh, aux, 11, 1, 0.3, True)
        ax.axis(axis)
        ax.set_title('overlap cells on task %d' % ir)
        fig.savefig(os.path.join(output_dir,
                                 'petsc_overlaps_%02d.png' % ir),
                    bbox_inches='tight')
        fig.delaxes(ax)

def plot_local_dofs(axs, field, field_i, omega_gi, output_dir, rank):
    """
    Plot the local ang global field DOFs local to the subdomain on the task
    with the given `rank`.
    """
    mesh = field.domain.mesh

    ax = pc.plot_wireframe(axs[0], mesh.cmesh)
    coors = field_i.get_coor()
    econn = field_i.econn
    ax = pd.plot_global_dofs(ax, coors, econn)
    ax.set_title('local DOFs on task %d' % rank)
    ax.figure.savefig(os.path.join(output_dir, 'local_dofs_%02d.png' % rank),
                      bbox_inches='tight')

    ax = pc.plot_wireframe(axs[1], mesh.cmesh)
    coors = field.get_coor()
    econn = field.get_econn(('cell', omega_gi.tdim), omega_gi, 0)
    ax = pd.plot_global_dofs(ax, coors, econn)
    ax.set_title('global DOFs on task %d' % rank)
    ax.figure.savefig(os.path.join(output_dir, 'local_global_%02d.png' % rank),
                      bbox_inches='tight')

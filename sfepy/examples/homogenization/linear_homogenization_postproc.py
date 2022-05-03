"""
This example shows how to use the VTK postprocessing functions.
"""

from __future__ import absolute_import
import os.path as osp
from .linear_homogenization import *
from sfepy.postprocess.utils_vtk import get_vtk_from_mesh,\
    get_vtk_by_group, get_vtk_surface, get_vtk_edges, write_vtk_to_file,\
    tetrahedralize_vtk_mesh

options.update({
    'post_process_hook' : 'post_process',
})

def post_process(out, problem, state, extend=False):

    mesh = problem.domain.mesh
    mesh_name = mesh.name[mesh.name.rfind(osp.sep) + 1:]

    vtkdata = get_vtk_from_mesh(mesh, out, 'postproc_')
    matrix = get_vtk_by_group(vtkdata, 1, 1)

    matrix_surf = get_vtk_surface(matrix)
    matrix_surf_tri = tetrahedralize_vtk_mesh(matrix_surf)
    write_vtk_to_file('%s_mat1_surface.vtk' % mesh_name, matrix_surf_tri)

    matrix_edges = get_vtk_edges(matrix)
    write_vtk_to_file('%s_mat1_edges.vtk' % mesh_name, matrix_edges)

    return out

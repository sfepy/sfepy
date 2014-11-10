"""Postprocessing utils based on VTK library"""

import vtk
import os.path as osp

def get_vtk_from_file(filename):
    """
    Read VTK file.

    Parameters
    ----------
    filename : str
        Name of the VTK file.

    Returns
    -------
    vtkdata : VTK object
        Mesh, scalar, vector and tensor data.
    """
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.ReadAllTensorsOn()
    reader.Update()

    return reader.GetOutput()

def write_vtk_to_file(filename, vtkdata):
    """
    Write VTK file.

    Parameters
    ----------
    filename : str
        Name of the VTK file.

    vtkdata : VTK object
        Mesh, scalar, vector and tensor data.
    """
    writer = vtk.vtkGenericDataObjectWriter()
    writer.SetFileName(filename)
    writer.SetInput(vtkdata)
    writer.Update()

def get_vtk_from_mesh(mesh, data, prefix=''):
    mesh_name = mesh.name[mesh.name.rfind(osp.sep) + 1:]
    vtkname = '%s%s.vtk' % (prefix, mesh_name)
    mesh.write(vtkname, io='auto', out=data)
    vtkdata = get_vtk_from_file(vtkname)

    return vtkdata

def get_vtk_surface(vtkdata):
    """
    Get mesh surface.

    Parameters
    ----------
    vtkdata : VTK object
        Mesh, scalar, vector and tensor data.

    Returns
    -------
    surface : VTK object
        Mesh, scalar, vector and tensor data.
    """
    surface = vtk.vtkDataSetSurfaceFilter()
    surface.SetInput(vtkdata)
    surface.Update()

    return surface.GetOutput()

def get_vtk_edges(vtkdata):
    """
    Get mesh edges.

    Parameters
    ----------
    vtkdata : VTK object
        Mesh, scalar, vector and tensor data.

    Returns
    -------
    edges : VTK object
        Mesh, scalar, vector and tensor data.
    """
    edges = vtk.vtkExtractEdges()
    edges.SetInput(vtkdata)
    edges.Update()

    return edges.GetOutput()

def get_vtk_by_group(vtkdata, group_lower, group_upper=None):
    """
    Get submesh by material group id.

    Parameters
    ----------
    vtkdata : VTK object
        Mesh, scalar, vector and tensor data.

    group_lower : int
        The lower material id.

    group_lower : int
        The Upper material id.

    Returns
    -------
    slection : VTK object
        Mesh, scalar, vector and tensor data.
    """
    selection = vtk.vtkThreshold()
    selection.SetInput(vtkdata)
    selection.SetInputArrayToProcess(0, 0, 0,
                                     vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,
                                     "mat_id")
    if group_upper is None:
        group_upper = group_lower

    selection.ThresholdBetween(group_lower, group_upper)
    selection.Update()

    return selection.GetOutput()

def tetrahedralize_vtk_mesh(vtkdata):
    """
    3D cells are converted to tetrahedral meshes, 2D cells to triangles.

    Parameters
    ----------
    vtkdata : VTK object
        Mesh, scalar, vector and tensor data.

    Returns
    -------
    tetra : VTK object
        Mesh, scalar, vector and tensor data.
    """
    tetra = vtk.vtkDataSetTriangleFilter()
    tetra.SetInput(vtkdata)
    tetra.Update()

    return tetra.GetOutput()

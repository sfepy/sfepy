#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Mesh Preview Generator.

Examples
--------

$ ./script/gen_mesh_prev.py meshes/2d/
"""
from __future__ import print_function
from __future__ import absolute_import
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import sys
sys.path.append('.')
import os
import vtk
from sfepy.discrete.fem import Mesh

def gen_shot(vtk_filename, png_filename):
    """
    Generate PNG image of the FE mesh.

    Parameters
    ----------
    vtk_filename : str
        The input mesh filename (file in VTK format).

    png_filename : str
        The name of the output PNG file.
    """

    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(vtk_filename)
    reader.Update()
    bnd = reader.GetOutput().GetPoints().GetBounds()

    surface0 = vtk.vtkDataSetSurfaceFilter()
    surface0.SetInput(reader.GetOutput())
    surface0.Update()

    if abs(bnd[5] - bnd[4]) > 1.0e-12:
        tr = vtk.vtkTransform()
        tr.RotateWXYZ(45,1,1,1)

        trFilter = vtk.vtkTransformPolyDataFilter()
        trFilter.SetTransform(tr)
        trFilter.SetInputConnection(surface0.GetOutputPort())
        trFilter.Update()
        surface = trFilter

    else:
        surface = surface0

    ca,cb = surface.GetOutput().GetCellData().GetScalars().GetRange()

    lut = vtk.vtkLookupTable()
    lut.SetHueRange(0.667, 0.667)
    lut.SetSaturationRange(0.0, 1.0)
    lut.SetValueRange(0.8, 1.0)
    lut.SetAlphaRange(1.0, 1.0)
    lut.SetTableRange(ca,cb)

    gf = vtk.vtkGraphicsFactory()
    gf.SetOffScreenOnlyMode(1)
    gf.SetUseMesaClasses(1)

    ifa = vtk.vtkImagingFactory()
    ifa.SetUseMesaClasses(1)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetLookupTable(lut)
    mapper.SetScalarRange(ca,cb);
    mapper.SetInput(surface.GetOutput())
    mapper.SetScalarModeToUseCellData()

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    mapper2 = vtk.vtkPolyDataMapper()
    mapper2.SetInput(surface.GetOutput())
    actor2 = vtk.vtkActor()
    actor2.SetMapper(mapper2)
    actor2.GetProperty().SetRepresentationToWireframe()

    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.SetOffScreenRendering(1)
    renWin.AddRenderer(ren)
    ren.AddActor(actor)
    ren.AddActor(actor2)
    renWin.Render()

    image = vtk.vtkWindowToImageFilter()
    image.SetInput(renWin)
    image.Update()

    base, _ = os.path.splitext(vtk_filename)
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(png_filename)
    writer.SetInput(image.GetOutput())
    writer.Write()

def main():
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('mesh_dir')
    options = parser.parse_args()

    mesh_dir = options.mesh_dir

    mesh_files = []
    for (dirpath, dirnames, filenames) in os.walk(mesh_dir):
        for ii in filenames:
            _, ext = os.path.splitext(ii)
            if ext.lower() in ['.mesh', '.vtk']:
                mesh_files.append(dirpath + os.path.sep + ii)

    for ii in mesh_files:
        base, ext = os.path.splitext(ii)
        fname_out = base + '.png'
        if ext == '.mesh':
            fname_in = 'aux.vtk'
            mesh = Mesh.from_file(ii)
            mesh.write(fname_in, io='auto')

        else:
            fname_in = ii

        print(('writing %s...' % fname_out))
        gen_shot(fname_in, fname_out)

if __name__ == "__main__":
    main()

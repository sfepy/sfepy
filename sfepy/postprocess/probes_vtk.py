"""Classes for probing values of Variables, for example, along a line,
using PyVTK library"""

from __future__ import absolute_import
import numpy as nm
import vtk
from vtk.util import numpy_support as vtknm
import os.path as osp

from sfepy.base.base import Struct, output
from sfepy.linalg import make_axis_rotation_matrix
from sfepy.postprocess.utils_vtk import get_vtk_from_mesh, get_vtk_from_file
from six.moves import range

vtk_version = vtk.vtkVersion().GetVTKMajorVersion()

class Probe(Struct):
    """
    Probe class.
    """

    def __init__(self, data, mesh, **kwargs):
        """
        Parameters
        ----------
        data : dict
            The output dictionary.
        mesh : Mesh
            The mesh.
        """

        Struct.__init__(self, name=mesh.name, **kwargs)

        self.mesh_name = mesh.name[mesh.name.rfind(osp.sep) + 1:]
        self.dim = mesh.dim
        self.vtkdata = get_vtk_from_mesh(mesh, data, 'probe_')

        self.vtkprobe = vtk.vtkProbeFilter()
        if vtk_version < 6:
            self.vtkprobe.SetSource(self.vtkdata)
        else:
            self.vtkprobe.SetSourceData(self.vtkdata)

        self.probes = {}
        self.probes_png = {}

    def new_vtk_polyline(self, points, closed=False):
        """
        Create the VTKPolyData object and store the line data.

        Parameters
        ----------
        points : array
            The line points.

        Returns
        -------
        vtkpd : VTK object
            VTKPolyData with the polyline.
        """

        npts = points.shape[0]
        pts = vtk.vtkPoints()
        pts.SetNumberOfPoints(npts)
        for ii in range(npts):
            pts.SetPoint(ii, points[ii,:])

        nlns = npts
        if closed:
            nlns += 1
        lns = vtk.vtkCellArray()
        lns.InsertNextCell(nlns)
        for ii in range(npts):
            lns.InsertCellPoint(ii)

        if closed:
            lns.InsertCellPoint(0)

        vtkpd = vtk.vtkPolyData()
        vtkpd.SetPoints(pts)
        vtkpd.SetLines(lns)
        if vtk_version < 6:
            vtkpd.Update()

        return vtkpd

    def add_line_probe(self, name, p0, p1, n_point):
        """
        Create the line probe - VTK object.

        Parameters
        ----------
        name : str
            The probe name.
        p0 : array_like
            The coordinates of the start point.
        p1 : array_like
            The coordinates of the end point.
        n_point : int
           The number of probe points.
        """

        line = vtk.vtkLineSource()
        line.SetPoint1(p0)
        line.SetPoint2(p1)
        line.SetResolution(n_point)
        line.Update()

        pars = nm.arange(n_point + 1) / float(n_point)
        self.probes[name] = (line, pars)
        self.probes_png[name] = False

    def add_points_probe(self, name, coors):
        """
        Create the point probe - VTK object.

        Parameters
        ----------
        name : str
            The probe name.
        coors : array
            The coordinates of the probe points.
        """
        coors = nm.asarray(coors)
        npts = coors.shape[0]
        pts = vtk.vtkPoints()
        pts.SetNumberOfPoints(npts)
        for ii in range(npts):
            pts.SetPoint(ii, coors[ii, :])

        poly = vtk.vtkPolyData()
        poly.SetPoints(pts)
        if vtk_version < 6:
            vtk.Update()

        pars = nm.arange(npts + 1)
        self.probes[name] = (poly, pars)
        self.probes_png[name] = False

    def add_ray_probe(self, name, p0, dirvec, p_fun, n_point):
        """
        Create the ray (line) probe - VTK object.

        Parameters
        ----------
        name : str
            The probe name.
        p0 : array
            The coordinates of the start point.
        dirvec : array
            The probe direction vector.
        p_fun : function
            The function returning the probe parametrization along the dirvec
            direction.
        n_point : int
           The number of probe points.
        """

        p0 = nm.array(p0, dtype=nm.float64)
        dirvec = nm.array(dirvec, dtype=nm.float64)
        dirvec /= nm.linalg.norm(dirvec)

        pars = p_fun(nm.arange(n_point, dtype=nm.float64))
        points = p0 + dirvec * pars[:,None]

        ray = self.new_vtk_polyline(points)
        self.probes[name] = (ray, pars)
        self.probes_png[name] = False

    def add_circle_probe(self, name, centre, normal, radius, n_point):
        """
        Create the ray (line) probe - VTK object.

        Parameters
        ----------
        name : str
            The probe name.
        centre : array
            The coordinates of the circle center point.
        normal : array
             The normal vector perpendicular to the circle plane.
        radius : float
            The radius of the circle.
        n_point : int
           The number of probe points.
        """

        pars = nm.linspace(0.0, 2.0*nm.pi, n_point + 1)[:-1]

        # Create the points in xy plane, centered at the origin.
        x = radius * nm.cos(pars[:,None])
        y = radius * nm.sin(pars[:,None])

        if len(centre) == 3:
            z = nm.zeros((n_point, 1), dtype=nm.float64)
            points = nm.c_[x, y, z]

            # Rotate to satisfy the normal, shift to the centre.
            n1 = nm.array([0.0, 0.0, 1.0], dtype=nm.float64)
            axis = nm.cross(n1, normal)
            angle = nm.arccos(nm.dot(n1, normal))

            if nm.linalg.norm(axis) < 0.1:
                # n1 == self.normal
                rot_mtx = nm.eye(3, dtype=nm.float64)
            else:
                rot_mtx = make_axis_rotation_matrix(axis, angle)

            points = nm.dot(points, rot_mtx)

        else:
            points = nm.c_[x, y]

        points += centre

        circle = self.new_vtk_polyline(points, closed=True)
        self.probes[name] = (circle, pars)
        self.probes_png[name] = False

    def gen_mesh_probe_png(self, probe, png_filename):
        """
        Generate PNG image of the FE mesh.

        Parameters
        ----------
        probe : VTK objectstr
            The probe, VTKPolyData or VTKSource.
        png_filename : str
            The name of the output PNG file.
        """

        surface = vtk.vtkDataSetSurfaceFilter()
        if vtk_version < 6:
            surface.SetInput(self.vtkdata)

        else:
            surface.SetInputData(self.vtkdata)

        surface.Update()

        gf = vtk.vtkGraphicsFactory()
        gf.SetOffScreenOnlyMode(1)
        gf.SetUseMesaClasses(1)

        ifa = vtk.vtkImagingFactory()
        ifa.SetUseMesaClasses(1)

        mapper = vtk.vtkPolyDataMapper()
        if vtk_version < 6:
            mapper.SetInput(surface.GetOutput())

        else:
            mapper.SetInputData(surface.GetOutput())

        mapper.SetScalarModeToUseCellData()

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetOpacity(0.33)

        mapper2 = vtk.vtkPolyDataMapper()
        if hasattr(probe, 'GetOutput'):
            probe0 = probe.GetOutput()

        else:
            probe0 = probe

        if vtk_version < 6:
            mapper2.SetInput(probe0)

        else:
            mapper2.SetInputData(probe0)

        actor2 = vtk.vtkActor()
        actor2.SetMapper(mapper2)
        actor2.GetProperty().SetColor(0,0,0)
        actor2.GetProperty().SetLineWidth(2)
        ren = vtk.vtkRenderer()
        renWin = vtk.vtkRenderWindow()
        renWin.SetOffScreenRendering(1)
        renWin.AddRenderer(ren)
        ren.AddActor(actor)
        ren.AddActor(actor2)
        ren.SetBackground(1, 1, 1)
        renWin.Render()

        image = vtk.vtkWindowToImageFilter()
        image.SetInput(renWin)
        image.Update()

        writer = vtk.vtkPNGWriter()
        writer.SetFileName(png_filename)
        if vtk_version < 6:
            writer.SetInput(image.GetOutput())

        else:
            writer.SetInputData(image.GetOutput())

        writer.Write()

    def __call__(self, probe_name, variable, probe_view=False,
                 ret_points=False):
        """
        Do the probe for the given variable.

        Parameters
        ----------
        probe_name : str
            The name of previously defined probe.
        variable : str
            The variable to be probed.
        probe_view: bool
            If True, save the probe visualization.
        ret_points : bool
            If True, return also the probe points.

        Returns
        -------
        params : array
            The parametrization of the probe points.
        points : array, optional
            If `ret_points` is True, the coordinates of points corresponding to
            `params`, where the `variable` is evaluated.
        values : array
            The probe values in the points.
        """

        inp = self.probes[probe_name][0]
        params = self.probes[probe_name][1]

        if hasattr(inp, 'GetOutputPort'):
            self.vtkprobe.SetInputConnection(inp.GetOutputPort())

        else:
            if vtk_version < 6:
                self.vtkprobe.SetInput(inp)

            else:
                self.vtkprobe.SetInputData(inp)

        self.vtkprobe.Update()
        pdata = self.vtkprobe.GetOutput()
        values = vtknm.vtk_to_numpy(pdata.GetPointData().GetArray(variable))

        output('probing data')
        output('  probe name: %s' % probe_name)
        output('  variable: %s' % variable)
        output('  points: %s' % params.shape[0])

        if probe_view and not(self.probes_png[probe_name]):
            pngname = 'probe_%s_%s.png' % (self.mesh_name, probe_name)
            self.gen_mesh_probe_png(inp, pngname)
            self.probes_png[probe_name] = True

        if ret_points:
            points = vtknm.vtk_to_numpy(pdata.GetPoints().GetData())
            points = points[:, :self.dim]

            return params, points, values

        else:
            return params, values


class ProbeFromFile(Probe):
    """
    Probe class - read a given VTK file.
    """

    def __init__(self, filename, **kwargs):
        """
        Parameters
        ----------
        filename : dict
            The name of a VTK file.
        """

        bname = osp.splitext(osp.basename(filename))[0]
        Struct.__init__(self, name=bname, **kwargs)

        self.vtkdata = get_vtk_from_file(filename)
        self.mesh_name = bname
        self.dim = 3

        self.vtkprobe = vtk.vtkProbeFilter()
        if vtk_version < 6:
            self.vtkprobe.SetSource(self.vtkdata)
        else:
            self.vtkprobe.SetSourceData(self.vtkdata)

        self.probes = {}
        self.probes_png = {}

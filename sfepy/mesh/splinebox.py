import numpy as nm
from bspline import BSpline

class SplineBox:
    """
    B-spline geometry parametrization. Geometry can be modified
    by moving spline control points.
    """
    @staticmethod
    def create_spb(bbox, coors, degree=3, nsg=None):
        bbox = nm.array(bbox)
        coors = nm.array(coors)
        dim = coors.shape[1]
        spbox = []

        if nsg is None:
            nsg = nm.ones((dim,), dtype=nm.int)
        else:
            nsg = nm.array(nsg)

        ncpoints = 1
        spbox = {'base': [None] * dim, 'uidx': [None] * dim,
                 'cp': [None] * dim, 'ncp': [None] * dim,
                 'cp_idx': [None] * dim}
        for idim in range(dim):
            ucoors, ucoors_idx = nm.unique(coors[:,idim], return_inverse=True)
            ncp = degree + nsg[idim]
            bspl = BSpline(degree, ncp=ncp)
            bspl.make_knot_vector(knot_range=(bbox[idim,0], bbox[idim,1]))
            knots = bspl.get_knot_vector()
            cp = nm.zeros((ncp,), dtype=nm.double)
            for j in range(cp.shape[0]):
                cp[j] = nm.sum(knots[(j + 1):(j + degree + 1)]) / degree

            spbox['base'][idim] = bspl.eval_basis(t=ucoors,
                                                  return_val=True).copy()
            spbox['uidx'][idim] = ucoors_idx[:]
            spbox['cp'][idim] = cp[:]
            spbox['ncp'][idim] = ncp
            ncpoints *= ncp

        cpoints = nm.zeros((ncpoints, dim), dtype=nm.double)
        ncp = spbox['ncp']
        cp = spbox['cp']
        if dim == 2:
            idxs = nm.mgrid[0:ncp[0],0:ncp[1]]

        elif dim == 3:
            idxs = nm.mgrid[0:ncp[0],0:ncp[1],0:ncp[2]]

        aux = [1]
        for ii in range(dim - 1):
            aux.append(aux[ii] * ncp[ii])
        spbox['mul_cp_idx'] = nm.array(aux)

        for ii in range(dim):
            spbox['cp_idx'][ii] = idxs[ii].reshape(ncpoints, order='F')
            cpoints[:,ii] = cp[ii][spbox['cp_idx'][ii]]

        return spbox, cpoints

    def __init__(self, bbox, coors,
                 name='spbox', **kwargs):
        """
        Create a SplineBox.

        Parameters
        ----------
        bbox : array
            Mesh bounding box.
        coors : array
            Coordinates of mesh nodes.
        name : str
            Object name.
        """
        self.spbox, self.control_points = self.create_spb(bbox, coors)
        self.ncoors, self.dim = coors.shape
        self.control_points0 = self.control_points.copy()

    def get_coors_shape(self):
        """
        Get the shape of the coordinates.
        """
        return (self.ncoors, self.dim)

    def get_control_points(self, init=False):
        """
        Get spline control points coordinates.

        Returns
        -------
        cpt_coors : array
            The coordinates of the spline control points.
        init : bool
            If True, return initial state.
        """
        if init:
            return self.control_points0

        else:
            return self.control_points

    def set_control_points(self, cpt_coors, add=False):
        """
        Set spline control points position.

        Parameters
        ----------
        cpt_coors : array
            The coordinates of the spline control points.
        add : bool
            If True, coors += cpt_coors
        """
        if add:
            self.control_points += cpt_coors

        else:
            self.control_points = cpt_coors

    def change_shape(self, cpoint, val):
        """
        Change shape of spline parametrization.

        Parameters
        ----------
        cpoint : list
            The indices of the spline control point.
        val : array
            Displacement.
        """
        idx = nm.dot(nm.array(cpoint), self.spbox['mul_cp_idx'])
        self.control_points[idx,:] += val

    def evaluate(self, cp_coors=None):
        """
        Evaluate SplineBox.

        Returns
        -------
        coors : array
            The coordinates corresponding to the actual spline control points
            position.
        cp_coors : array
            If is not None, use as control points cooardinates.
        """
        if cp_coors is None:
            cp_coors = self.control_points

        ncp, dim = cp_coors.shape
        base = self.spbox['base']
        uidx = self.spbox['uidx']
        cp_idx = self.spbox['cp_idx']
        ncoors = uidx[0].shape[0]
        coors = nm.zeros((ncoors, dim), dtype=nm.double)
        aux = nm.ones((ncoors,ncp), dtype=nm.double)
        for ii in range(self.dim):
            aux *= base[ii][uidx[ii],:][:,cp_idx[ii]]

        coors = nm.dot(aux, cp_coors)

        return coors

    def dvelocity(self, cpoint, dirvec):
        """
        Evaluate derivative of spline in a given control point and direction.

        Parameters
        ----------
        cpoint : list
            The indices of the spline control point.
        dir : array
            The directional vector.

        Returns
        -------
        dvel : array
            The design velocity field.
        """
        base = self.spbox['base']
        uidx = self.spbox['uidx']

        aux = nm.ones((self.ncoors,), dtype=nm.double)
        for ii in range(self.dim):
            aux *= base[ii][uidx[ii],cpoint[ii]]

        dvel = nm.dot(aux[:,nm.newaxis], nm.reshape(dirvec, (1,self.dim)))

        return dvel

    def write_vtk(self, filename):
        """
        Write the SplineBox shape to the VTK file.

        Parameters
        ----------
        filename : str
            The VTK file name.
        """
        ncp = self.spbox['ncp']
        npt = nm.prod(ncp)

        f = open(filename, 'w')
        f.write("# vtk DataFile Version 2.6\nspbox file\n"
                "ASCII\nDATASET UNSTRUCTURED_GRID\n\n")
        f.write("POINTS %d float\n" % npt)

        if self.dim == 2:
            ptformat = "%e %e 0.0\n"

        elif self.dim == 3:
            ptformat = "%e %e %e\n"

        for cpt in self.control_points:
            f.write(ptformat % tuple(cpt))

        cells = nm.array([nm.arange(0, ncp[0] - 1), nm.arange(1, ncp[0])]).T
        cp = ncp[0]
        nc = cp - 1
        for ii in range(1, self.dim):
            cells1 = []
            ncpi = ncp[ii]
            for jj in range(ncpi):
                cells1.append(cells + jj * cp)
            nc = nc * ncpi

            cells = nm.array([nm.arange(0, ncpi - 1),
                              nm.arange(1, ncpi)]).T * cp
            for jj in range(cp):
                cells1.append(cells + jj)
            nc += (ncpi - 1) * cp

            cells = nm.vstack(cells1)
            cp *= ncp[ii]

        f.write("\nCELLS %d %d\n" % (nc, 3 * nc))
        for ii in cells:
            f.write("2 %d %d\n" % tuple(ii))
        f.write("\nCELL_TYPES %d\n" % nc)
        f.write("3\n" * nc)
        f.close()

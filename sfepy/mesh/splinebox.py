import numpy as nm
from bspline import BSpline
from sfepy.base.base import Struct

class SplineBox(Struct):
    """
    B-spline geometry parametrization. The geometry can be modified
    by moving spline control points.
    """
    @staticmethod
    def gen_cp_idxs(ncp):
        dim = len(ncp)

        if dim == 2:
            idxs = nm.mgrid[0:ncp[0],0:ncp[1]]
        elif dim == 3:
            idxs = nm.mgrid[0:ncp[0],0:ncp[1],0:ncp[2]]
        else:
            raise(ValueError)

        cp_idxs = []
        mul_idx = [1]
        for ii in range(dim - 1):
            mul_idx.append(mul_idx[ii] * ncp[ii])

        cp_idxs = []
        for ii in range(dim):
            cp_idxs.append(idxs[ii].reshape(nm.prod(ncp), order='F'))

        return cp_idxs, nm.array(mul_idx)

    @staticmethod
    def create_spb(bbox, coors, degree=3, nsg=None):
        nc, dim = coors.shape
        inside = nm.ones((nc,), dtype=nm.bool)
        nsg = nm.ones((dim,), dtype=nm.int) if nsg is None else nm.array(nsg)

        for idim in range(dim):
            inrange = nm.logical_and(coors[:,idim] >= bbox[idim][0],
                                     coors[:,idim] <= bbox[idim][1])
            inside = nm.logical_and(inside, inrange)

        ncpoints = 1
        base, uidx, ncp, cp = [], [], [], []
        for idim in range(dim):
            ucoors, ucoors_idx = nm.unique(coors[inside,idim],
                                           return_inverse=True)
            ncp0 = degree + nsg[idim]
            bspl = BSpline(degree, ncp=ncp0)
            bspl.make_knot_vector(knot_range=(bbox[idim][0], bbox[idim][1]))
            knots = bspl.get_knot_vector()
            cp0 = nm.zeros((ncp0,), dtype=nm.double)
            for j in range(cp0.shape[0]):
                cp0[j] = nm.sum(knots[(j + 1):(j + degree + 1)]) / degree

            base.append(bspl.eval_basis(t=ucoors, return_val=True))
            uidx.append(ucoors_idx)
            ncp.append(ncp0)
            cp.append(cp0)
            ncpoints *= ncp0

        cpoints = nm.zeros((ncpoints, dim), dtype=nm.double)
        cp_idx, mul_cp_idx = SplineBox.gen_cp_idxs(ncp)
        for ii in range(dim):
            cpoints[:,ii] = cp[ii][cp_idx[ii]]

        return {'base': base,
                'uidx': uidx,
                'ncp': ncp,
                'cp_idx': cp_idx,
                'mul_cp_idx': mul_cp_idx,
                'control_points': cpoints,
                'idxs_inside' : inside}

    def __init__(self, bbox, coors, nsg=None):
        """
        Create a SplineBox.

        Parameters
        ----------
        bbox : array
            The mesh bounding box.
        coors : array
            The coordinates of mesh nodes.
        nsg : array
            The number of segments.
        """
        bbox = nm.array(bbox)
        coors = nm.array(coors)
        self.__dict__.update(self.create_spb(bbox, coors, nsg=nsg))
        self.ncoors, self.dim = coors.shape
        self.coors = coors.copy()
        self.control_points0 = self.control_points.copy()

    def get_coors_shape(self):
        """
        Get the shape of the coordinates.
        """
        return (self.ncoors, self.dim)

    def get_control_points(self, init=False):
        """
        Get the spline control points coordinates.

        Returns
        -------
        cpt_coors : array
            The coordinates of the spline control points.
        init : bool
            If True, return the initial state.
        """
        if init:
            return self.control_points0

        else:
            return self.control_points

    def set_control_points(self, cpt_coors, add=False):
        """
        Set the spline control points position.

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

    def move_control_point(self, cpoint, val):
        """
        Change shape of spline parametrization.

        Parameters
        ----------
        cpoint : int, list
            The position (index or grid indicies) of the spline control point.
        val : array
            Displacement.
        """
        if type(cpoint) in [list, tuple, nm.ndarray]:
            idx = nm.dot(nm.array(cpoint), self.mul_cp_idx)
        else:
            idx = cpoint
        self.control_points[idx,:] += val

    def evaluate(self, cp_coors=None, outside=True):
        """
        Evaluate the new position of the mesh coordinates.

        Parameters
        ----------
        coors : array
            The coordinates corresponding to the actual spline control point
            positions.
        cp_coors : array
            If is not None, use as control point coordinates.
        outside : bool
            If True, return also the coordinates outside the spline box.

        Returns
        -------
        new_coors : array
            The new position of the mesh coordinates.
        """
        if cp_coors is None:
            cp_coors = self.control_points

        ncp, dim = cp_coors.shape
        aux = nm.ones((self.uidx[0].shape[0], ncp), dtype=nm.double)
        for ii in range(self.dim):
            aux *= self.base[ii][self.uidx[ii],:][:,self.cp_idx[ii]]

        if outside and hasattr(self, 'idxs_inside'):
            coors = self.coors.copy()
            coors[self.idxs_inside,:] = nm.dot(aux, cp_coors)
            return coors
        else:
            return nm.dot(aux, cp_coors)

    def evaluate_derivative(self, cpoint, dirvec):
        """
        Evaluate derivative of the spline
        in a given control point and direction.

        Parameters
        ----------
        cpoint : int, list
            The position (index or grid indicies) of the spline control point.
        dirvec : array
            The directional vector.

        Returns
        -------
        diff : array
            The derivative field.
        """
        if type(cpoint) in [list, tuple, nm.ndarray]:
            idxs = cpoint
        else:
            idxs = []
            aux = cpoint
            for ii in range(self.dim):
                idxs.append(aux / self.mul_cp_idx[-ii])
                aux = aux % self.mul_cp_idx[-ii]
            idxs = idxs[::-1]

        aux = nm.ones((self.uidx[0].shape[0],), dtype=nm.double)
        for ii in range(self.dim):
            aux *= self.base[ii][self.uidx[ii],cpoint[ii]]

        return nm.dot(aux[:,nm.newaxis], nm.reshape(dirvec, (1, self.dim)))

    def write_control_net(self, filename):
        """
        Write the SplineBox shape to the VTK file.

        Parameters
        ----------
        filename : str
            The VTK file name.
        """
        ncp = self.ncp
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

class SplineRegion2D(SplineBox):
    """
    B-spline geometry parametrization. The boundary of the SplineRegion2D
    is defined by BSpline curves.
    """
    @staticmethod
    def points_in_poly(points, poly, tol=1e-6):
        """
        Find which points are located inside the polygon.
        """
        poly = nm.array(poly)
        points = nm.array(points)
        inside = nm.zeros((points.shape[0],), dtype=nm.bool)

        p1 = poly[:-1]
        p2 = poly[1:]
        a1 = (p2[:,1] - p1[:,1])
        a1nz = nm.where(nm.fabs(a1) > 1e-16)[0]
        a2 = (p2[a1nz,0] - p1[a1nz,0]) / a1[a1nz]
        for jj, pt in enumerate(points):
            # on edges?
            if nm.any(nm.linalg.norm(p1 - pt, axis=1)
                      + nm.linalg.norm(p2 - pt, axis=1)
                      - nm.linalg.norm(p1 - p2, axis=1) < tol):
                inside[jj] = True
                continue

            # inside?
            val = nm.logical_and((p1[a1nz,1] > pt[1]) != (p2[a1nz,1] > pt[1]),
                                 pt[0] < (a2*(pt[1] - p1[a1nz,1]) + p1[a1nz,0]))

            if (nm.where(val)[0].shape[0] % 2) > 0:
                inside[jj] = True

        return nm.where(inside)[0]

    @staticmethod
    def define_control_points(cp_bnd_coors, ncp):
        """
        Find positions of "inner" control points depending on boundary splines.
        """
        nx, ny = ncp
        grid = nm.zeros(ncp, dtype=nm.int32)
        grid.T.flat = nm.arange(nx * ny)

        coors = nm.zeros((nx * ny, 2), dtype=nm.float64)
        idxs1 = nm.arange(nx)
        idxs2 = nm.arange(1, ny - 1) * nx
        bcnd = nm.hstack([idxs1, idxs2 + nx - 1,
                          idxs1[::-1] + nx * (ny - 1), idxs2[::-1]])
        coors[bcnd,:] = cp_bnd_coors

        for ii in range(1, nx - 1):
            for jj, t in enumerate(nm.linspace(0, 1, ny)[1:-1]):
                c = (1 - t) * coors[ii,:] + t * coors[ii + (ny - 1)*nx,:]
                coors[ii + nx*(jj + 1),:] = c

        inside = grid[1:-1,1:-1].flatten()
        for iiter in range(5):
            for ii in inside:
                dx = nm.array([0., 0.])
                for jj in [-1, +1, -nx, + nx]:
                    dx -= 0.25 * (coors[ii,:] - coors[ii + jj,:])
                coors[ii] += 0.1 * dx

        return coors

    @staticmethod
    def create_spb(spl_bnd, coors, rho=10):
        """
        Initialize SplineBox knots, control points, base functions, ...
        """
        dim = 2
        if coors.shape[1] != dim:
            print 'Only 2D SplineBoxSBnd is supported!'
            raise(ValueError)

        bnd_poly = []
        bnd_cp = []
        for s in spl_bnd:
            s.set_param_n(rho)
            bnd_poly.append(s.eval()[:-1])
            bnd_cp.append(s.get_control_points()[:-1,:])

        bnd_poly.append(bnd_poly[0][0,:])
        ncpoints = 1
        base, bspl, uidx, ncp =  [], [], [], []
        for idim, si in enumerate([0, 1]):
            s = spl_bnd[si]
            bspl0 = BSpline(s.degree, ncp=s.ncp)
            bspl0.set_knot_vector(s.knots)
            bspl.append(bspl0)
            base.append(None)
            uidx.append(None)
            ncp.append(s.ncp)
            ncpoints *= s.ncp

        cp_idx, mul_cp_idx = SplineBox.gen_cp_idxs(ncp)
        cpoints = SplineRegion2D.define_control_points(nm.vstack(bnd_cp), ncp)
        idxs_inside = SplineRegion2D.points_in_poly(coors, nm.vstack(bnd_poly))

        return {'base': base,
                'bspl': bspl,
                'uidx': uidx,
                'ncp': ncp,
                'cp_idx': cp_idx,
                'mul_cp_idx': mul_cp_idx,
                'control_points': cpoints,
                'idxs_inside': idxs_inside}

    def find_ts(self, coors):
        """
        Function finds parameters (t, s) corresponding to given points (coors).
        """
        from scipy.optimize import minimize

        def ptdist(x, coors, spb):
            for ii in range(spb.dim):
                spb.base[ii] = spb.bspl[ii].eval_basis(t=x[ii], return_val=True)

            coors_approx = spb.evaluate(outside=False)
            return nm.linalg.norm(coors - coors_approx)

        def gen_grid(spb, rho):
            grid = nm.mgrid[0:rho,0:rho]
            t = nm.linspace(0, 1, rho)
            for ii in range(spb.dim):
                spb.uidx[ii] = grid[ii,:].reshape(rho**self.dim, order='F')
                spb.base[ii] = spb.bspl[ii].eval_basis(t=t, return_val=True)

            return spb.evaluate(outside=False)

        rho = 100
        grid = gen_grid(self, rho)

        for ii in range(self.dim):
            self.uidx[ii] = nm.array([0])

        ts = nm.zeros((coors.shape[0], self.dim), dtype=nm.float64)
        for ii, ic in enumerate(coors):
            idx = nm.argmin(nm.linalg.norm(grid - ic, axis=1))
            x0 = nm.array([idx % rho, idx / rho]) / (rho - 1.)
            fce = lambda x: ptdist(x, ic, self)
            ts[ii] = minimize(fce, x0, method='nelder-mead',
                              options={'xtol': 1e-5, 'disp': False}).x

        return ts

    def __init__(self, spl_bnd, coors, rho=1e3):
        """
        Create a SplineBox which boundary is defined by B-spline curves.

        Parameters
        ----------
        spl_bnd : list
            The list of BSpline objects (counterclockwise)
            defining the SplineBox boundary.
        coors : array
            The coordinates of the mesh nodes.
        rho : float
            The density of points defining the boundary polygon.
        """
        coors = nm.array(coors)
        self.__dict__.update(self.create_spb(spl_bnd, coors, rho))
        self.ncoors, self.dim = coors.shape
        self.coors = coors.copy()

        self.ts = self.find_ts(coors[self.idxs_inside,:])
        self.control_points0 = self.control_points.copy()

        for idim in range(self.dim):
            ucoors, ucoors_idx = nm.unique(self.ts[:,idim], return_inverse=True)
            self.base[idim] = self.bspl[idim].eval_basis(t=ucoors,
                                                         return_val=True)
            self.uidx[idim] = ucoors_idx

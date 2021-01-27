from __future__ import print_function
from __future__ import absolute_import
import sys
from six.moves import range
sys.path.append('.')

import numpy as nm
from sfepy.base.base import Struct

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

nm_f64_eps = nm.finfo(nm.float64).eps

def to_ndarray(a):
    if a is None:
        return None
    else:
        a = nm.asarray(a)
        if len(a.shape) == 0:
            a = a.reshape(1)
        return a

class BSpline(Struct):
    """
    B-spline curve representation
    """

    def __init__(self, degree=3, is_cyclic=False, ncp=0):
        """
        Initialize B-spline class.

        Parameters
        ----------
        degree : int
            The degree of the spline function.
        is_cyclic : bool
            Cyclic spline?.
        ncp : int
            The number of control points.
        """
        self.degree = degree
        self.knot_type = None
        self.is_cyclic = is_cyclic
        self.ncp = ncp
        self.knots = None
        self.basis = None
        self.curve_coors = None
        self.cp_coors = None
        self.approx_coors = None
        self.t = None

    def set_control_points(self, coors, cyclic_form=False):
        """
        Set the B-spline control points.

        Parameters
        ----------
        coors : array
            The coordinates of unique control points.
        cyclic_form : bool
            Are the control points in the cyclic form?
        """
        coors = to_ndarray(coors)

        if self.is_cyclic and not cyclic_form:
            coors = nm.vstack((coors, coors[0:self.degree,:]))

        self.cp_coors = coors
        self.ncp = coors.shape[0]

    def get_control_points(self):
        """
        Get the B-spline control points.

        Returns
        -------
        coors : array
            The coordinates of control points.
        """
        if self.is_cyclic:
            return self.cp_coors[:-self.degree,:]
        else:
            return self.cp_coors

    def set_param(self, t):
        """
        Set the B-spline parametric vector.

        Parameters
        ----------
        t : array
            The parameter vector of the B-spline.
        """
        self.t = to_ndarray(t)

        if self.knots is not None:
            endval = self.knots[-(self.degree + 1)]
            idxs = nm.where(self.t == endval)[0]
            self.t[idxs] -= nm_f64_eps

    def set_param_n(self, n=100, knot_range=(0.0, 1.0)):
        """
        Generate the B-spline parametric vector using the number of steps.

        Parameters
        ----------
        n : array
            The number of steps in the B-spline parametric vector.
        """
        self.t = nm.linspace(knot_range[0], knot_range[1], n)
        self.t[-1] -= nm_f64_eps

    @staticmethod
    def basis_function_dg0(t, knots, n):
        """
        Basis function: degree = 0

        Parameters
        ----------
        t : array
            The parametric vector.
        knots : array
            The knot vector.
        n : int
            The number of intervals.

        Returns
        -------
        bfun : array
           The spline basis function evaluated for given values.
        """
        nt = len(t)
        bfun = nm.zeros((nt,n), dtype=nm.float64)
        for ii in range(n):
            idxs = nm.where(nm.logical_and(knots[ii] <= t,
                                           t < knots[ii + 1]))[0]
            bfun[idxs,ii] = 1.0

        return bfun

    @staticmethod
    def basis_function_dg(degree, t, knots, n):
        """
        B-spline basis functions.

        Parameters
        ----------
        degree : int
            The degree of the spline function.
        t : array
            The parametric vector.
        knots : array
            The knot vector.
        n : int
            The number of intervals.

        Returns
        -------
        bfun : array
           The spline basis function evaluated for given values.
        """
        if degree >= 1:
            bfun_dgm1 = BSpline.basis_function_dg(degree - 1, t,
                                                  knots, n + 1)

            nt = len(t)
            bfun = nm.zeros((nt,n), dtype=nm.float64)
            for ii in range(n):
                c1 = t - knots[ii]
                c2 = knots[ii + degree] - knots[ii]

                if nm.abs(c2) > nm_f64_eps:
                     bfun[:,ii] = c1 / c2 * bfun_dgm1[:,ii]

                c1 = knots[ii + degree + 1] - t
                c2 = knots[ii + degree + 1] - knots[ii + 1]
                if nm.abs(c2) > nm_f64_eps:
                    bfun[:,ii] += c1 / c2 * bfun_dgm1[:,ii + 1]

        else:
            bfun = BSpline.basis_function_dg0(t, knots, n)

        return bfun

    def make_knot_vector(self, knot_type='clamped', knot_data=None,
                         knot_range=(0.0, 1.0)):
        """
        Create a knot vector of the requested type.

        Parameters
        ----------
        knot_type : str
            The knot vector type: clamped/cyclic/userdef.
        knot_data :
            The extra knot data.
        """
        if self.is_cyclic and 'cyclic' not in knot_type:
            knot_type = 'cyclic'

        ncp = self.ncp
        dg = self.degree
        n_knots = dg + ncp + 1
        n_inter = n_knots - 2 * dg
        aux = nm.linspace(knot_range[0], knot_range[1], n_inter)
        if knot_type == '' or knot_type == 'cyclic':
            dd = aux[1]
            self.knots = nm.hstack((nm.arange(-dg, 0) * dd,
                                    aux,
                                    nm.arange(1, dg + 1) * dd + 1))

        elif knot_type == 'clamped':
            self.knots = nm.array([aux[0]]* dg + list(aux) + [aux[-1]]* dg,
                                  dtype=nm.float64)

        else:
            raise NotImplementedError

        self.knot_type = knot_type

    def set_knot_vector(self, knots):
        """
        Set the knot vector.

        Parameters
        ----------
        knots : array
            The knot vector.
        """
        self.knot_type = 'userdef'
        self.knots = to_ndarray(knots)

    def get_knot_vector(self):
        """
        Return the knot vector.

        Returns
        -------
        knots : array
            The knot vector.
        """
        return self.knots

    def insert_knot(self, new):
        """
        Insert a new knot into the knot vector.

        Parameters
        ----------
        new : float
            The new knot value.
        """
        kn = self.knots
        dg = self.degree
        ncp = self.ncp
        cp = self.cp_coors
        idx = nm.where(nm.logical_and(kn[:-1] <= new, new < kn[1:]))[0]

        if len(idx) > 0:
            multi = len(nm.where(kn == new)[0])
            if multi < dg:
                # new knot
                newkn = nm.zeros((len(kn) + 1,), dtype=nm.float64)
                newkn[:(idx + 1)] = kn[:(idx + 1)]
                newkn[idx + 1] = new
                newkn[(idx + 2):] = kn[(idx + 1):]
                u1 = idx - dg + 1
                u2 = idx + 1

                # new control points
                newcp = nm.zeros((ncp + 1, cp.shape[1]), dtype=nm.float64)
                newcp[:u1,:] = cp[:u1,:]
                newcp[u2:,:] = cp[(u2 - 1):,:]

                for ii in range(u1, u2):
                    kn1 = kn[ii]
                    kn2 = kn[ii + dg]
                    dd = kn2 - kn1
                    newcp[ii,:] = (kn2 - new) / dd * cp[ii - 1] + \
                                  (new - kn1) / dd * cp[ii]

                self.knots = newkn
                self.cp_coors = newcp
                self.ncp = newcp.shape[0]

                # evaluate B-spline base functions for new configuration
                self.eval_basis()

            else:
                print('knot insertion failed: multiplicity = spline degree!')

        else:
            print('knot insertion failed: out of bounds!')

    def eval_basis(self, t=None, return_val=False):
        """
        Evaluate the basis of the bpsline.

        Parameters
        ----------
        t : array
            The parameter vector of the B-spline.
        """
        if t is not None:
            self.set_param(t)

        if self.knots is None:
            self.make_knot_vector()

        if self.t is None:
            self.set_param_n()

        self.basis = self.basis_function_dg(self.degree, self.t,
                                            self.knots, self.ncp)

        if return_val:
            return self.basis

    def eval(self, t=None, cp_coors=None):
        """
        Evaluate the coordinates of the bpsline curve.

        Parameters
        ----------
        t : array
            The parameter vector of the B-spline.
        cp_coors : array
            The coordinates of the control points.
        """
        if cp_coors is not None:
            self.set_control_points(cp_coors)

        self.eval_basis(t)
        self.curve_coors = nm.dot(self.basis, self.cp_coors)

        return self.curve_coors

    def draw(self, ret_ax=False, ax=None, color='r', cp_id=True):
        """
        Draw B-spline curve.

        Parameters
        ----------
        ret_ax : bool
            Return an axes object?
        ax : axes object
            The axes to which will be drawn.
        color : str
            Line color.
        cp_id : bool
            If True, label control points.
        """
        if self.curve_coors is None:
            self.eval()

        cc = self.curve_coors
        cp = self.cp_coors
        ci = self.approx_coors
        if ci is not None and self.is_cyclic:
            ci = nm.vstack((ci, ci[0,:]))

        if cc.shape[1] == 3:
            if ax is None:
                fig = plt.figure()
                ax = Axes3D(fig)
            ax.plot(cc[:,0], cc[:,1], cc[:,2], color + '-')
            if cp_id:
                ax.plot(cp[:,0], cp[:,1], cp[:,2], 'ko:', alpha=0.6)
            if ci is not None:
                ax.plot(ci[:,0], ci[:,1], ci[:,2], 'b--', alpha=0.6)

        else:
            if ax is None:
                fig = plt.figure()
                ax = fig.add_subplot(111)
            ax.plot(cc[:,0], cc[:,1], color + '-')
            if cp_id:
                ax.plot(cp[:,0], cp[:,1], 'ko:', alpha=0.6)
                for jj, icp in enumerate(self.cp_coors):
                    ax.text(icp[0], icp[1], 'N%d' % (jj + 1), fontsize=10)
            if ci is not None:
                ax.plot(ci[:,0], ci[:,1], 'b--', alpha=0.6)

        ax.set_aspect('equal')

        if ret_ax:
            return ax

        else:
            plt.show()

    def draw_basis(self):
        """
        Draw B-spline curve.
        """
        plt.figure()
        plt.plot(self.t,self.basis)
        plt.legend(['b%d' % (ii + 1) for ii in range(self.basis.shape[1])])
        plt.show()

    def approximate(self, coors, ncp=None, knot_type='clamped',
                    knots=None, alpha=0.5,
                    do_eval=False, do_param_correction=False):
        """
        Approximate set of points by the B-spline curve.

        Parameters
        ----------
        coors : array
            The coordinates of the approximated points.
        ncp : int
            The number of control points.
        knot_type : str
            The knot vector type.
        knots : array
            The knot vector.
        alpha : float
            The parameter vector distribution:
                1.0 = chordal
                0.5 = centripetal
        do_eval : bool
            Evaluate the curve coordinates?
        do_param_correction : bool
            Perform parametric corrections to improve the approximation?
        """
        coors = to_ndarray(coors)
        dg = self.degree

        if ncp is not None:
            if self.is_cyclic:
                ncp += dg
            self.ncp = ncp
            self.make_knot_vector(knot_type)

        if knots is not None:
            self.knots = knots
            self.knot_type = 'userdef'
            ncp = len(knots) - dg - 1
            self.ncp = ncp

        # param vector
        dist = nm.linalg.norm(coors[:-1,:] - coors[1:,:], axis=1)
        dista = nm.power(dist, alpha)
        self.t = nm.zeros((coors.shape[0],), dtype=nm.float64)
        self.t[1:] += dista.cumsum() / dista.sum()
        self.t[-1] -= nm_f64_eps

        while True:
            self.basis = self.basis_function_dg(dg, self.t,
                                                self.knots, ncp)

            A = nm.dot(self.basis.T, self.basis)
            b = nm.dot(self.basis.T, coors)
            # cyclic spline
            if self.is_cyclic:
                nred = ncp - dg
                R = nm.zeros((ncp, nred), dtype=nm.float64)
                for ii in range(nred):
                    R[ii,ii] = 1.0

                for ii in range(self.degree):
                    R[nred + ii,ii] = 1.0

                A = nm.dot(R.T, nm.dot(A, R))
                b = nm.dot(R.T, b)

                self.cp_coors = nm.dot(R, nm.dot(nm.linalg.inv(A), b))

            else:
                self.cp_coors = nm.dot(nm.linalg.inv(A), b)

            self.approx_coors = coors

            if not do_param_correction:
                break

        if do_eval:
            self.curve_coors = nm.dot(self.basis, self.cp_coors)

    def set_approx_points(self, coors):
        """
        Set the coordinates of approximated points.

        Parameters
        ----------
        coors : array
            The coordinates of approximated points.
        """
        self.approx_coors = to_ndarray(coors)

class BSplineSurf(Struct):
    """
    B-spline surface representation
    """

    def __init__(self, degree=(3,3), is_cyclic=(False, False)):
        """
        Initialize B-spline class.

        Parameters
        ----------
        degree : tuple of int
            The degree of the spline functions.
        is_cyclic : tuple of bool
            Cyclic splines?.
        """

        self.splines = [None, None]
        for ii in range(2):
            self.splines[ii] = BSpline(degree[ii], is_cyclic=is_cyclic[ii])

        self.surf_coors = None
        self.cp_coors = None
        self.approx_coors = None

    def set_control_points(self, coors, cyclic_form=False):
        """
        Set the B-spline control points.

        Parameters
        ----------
        coors : array
            The coordinates of unique control points.
        cyclic_form : bool
            Are the control points in the cyclic form?
        """
        coors = to_ndarray(coors)

        if self.splines[0].is_cyclic and not cyclic_form:
            coors = nm.vstack((coors, coors[0:self.splines[0].degree,:,:]))

        if self.splines[1].is_cyclic and not cyclic_form:
            coors = nm.hstack((coors, coors[:,0:self.splines[1].degree,:]))

        self.cp_coors = coors
        for ii in range(2):
            self.splines[ii].ncp = coors.shape[ii]

    def get_control_points(self):
        """
        Get the B-spline surface control points.

        Returns
        -------
        coors : array
            The coordinates of control points.
        """

        aux = self.cp_coors
        if self.splines[0].is_cyclic:
            aux = aux[:-self.splines[0].degree,:,:]

        if self.splines[1].is_cyclic:
            aux = aux[:,:-self.splines[1].degree,:]

        return aux

    def make_knot_vector(self, knot_type=('clamped', 'clamped'),
                         knot_data=(None, None)):
        """
        Create a knot vector of the requested type.

        Parameters
        ----------
        knot_type : tuple of str
            The knot vector types.
        knot_data : tuple of ANY
            The extra knot data.
        """
        for ii in range(2):
            self.splines[ii].make_knot_vector(knot_type[ii], knot_data[ii])

    def set_param_n(self, n=(100, 100)):
        """
        Generate the B-spline parametric vector using the number of steps.

        Parameters
        ----------
        n : tuple of array
            The number of steps in the B-spline parametric vectors.
        """
        for ii in range(2):
            self.splines[ii].set_param_n(n[ii])

    def set_approx_points(self, coors):
        """
        Set the coordinates of approximated points.

        Parameters
        ----------
        coors : array
            The coordinates of approximated points.
        """
        self.approx_coors = to_ndarray(coors)

    def eval(self, t=(None, None), cp_coors=None):
        """
        Evaluate the coordinates of the bpsline curve.

        Parameters
        ----------
        t : tuple of array
            The parametric vector of the B-splines.
        cp_coors : array
            The coordinates of the control points.
        """
        if cp_coors is not None:
            self.set_control_points(cp_coors)

        for ii in range(2):
            self.splines[ii].eval_basis(t[ii])

        nt = (len(self.splines[0].t), len(self.splines[1].t))
        ncp = (self.splines[0].ncp, self.splines[1].ncp)
        aux = nm.zeros((nt[0], ncp[1], 3), dtype=nm.float64)
        for ii in range(ncp[1]):
            aux[:,ii,:] = nm.dot(self.splines[0].basis, self.cp_coors[:,ii,:])

        self.surf_coors = nm.zeros(nt + (3,), dtype=nm.float64)
        for ii in range(nt[0]):
            self.surf_coors[ii,:,:] = nm.dot(self.splines[1].basis, aux[ii,:,:])

        return self.surf_coors

    def draw(self, ret_ax=False, ax=None):
        """
        Draw B-spline surface.

        Parameters
        ----------
        ret_ax : bool
            Return an axes object?
        ax : axes object
            The axes to which will be drawn.
        """

        if self.surf_coors is None:
            self.eval()

        fig = plt.figure()
        ax = Axes3D(fig)
        coors = self.surf_coors
        cs = coors.shape
        for ii in range(cs[0] - 1):
            for jj in range(cs[1] - 1):
                verts = nm.array([coors[ii,jj,:],
                                  coors[ii,jj + 1,:],
                                  coors[ii + 1,jj + 1,:],
                                  coors[ii + 1,jj,:]])

                quad = Poly3DCollection([verts],
                                        facecolors='gray', edgecolor='k',
                                        linewidths=0.2, alpha=0.5)
                ax.add_collection3d(quad)

        cp = self.cp_coors
        for ii in range(cp.shape[1]):
            ax.plot(cp[:,ii,0], cp[:,ii,1], cp[:,ii,2], 'ro--', linewidth=2.0)
        for ii in range(cp.shape[0]):
            ax.plot(cp[ii,:,0], cp[ii,:,1], cp[ii,:,2], 'ro--', linewidth=2.0)

        ax.set_aspect('equal')
        plt.show()

    def approximate(self, coors, ncp, do_eval=False):
        """
        Approximate set of points by the B-spline surface.

        Parameters
        ----------
        coors : array
            The coordinates of the approximated points.
        ncp : tuple of int
            The number of control points.
        """
        coors = to_ndarray(coors)

        nsc = coors.shape[0:2]
        aux = nm.zeros((nsc[0], ncp[1], 3), dtype=nm.float64)
        spl1 = self.splines[1]
        for ii in range(nsc[0]):
            spl1.approximate(coors[ii,...], ncp[1])
            aux[ii,...] = spl1.get_control_points()

        self.cp_coors = nm.zeros((ncp[0], ncp[1], 3), dtype=nm.float64)
        spl2 = self.splines[0]
        for ii in range(ncp[1]):
            spl2.approximate(aux[:,ii,:], ncp[0])
            self.cp_coors[:,ii,:] = spl2.get_control_points()

        self.approx_coors = coors

    def write_surface_vtk(self, filename, float_format='%.6f'):
        """
        Write the spline surface to VTK file.

        Parameters
        ----------
        filename: str
            Name of the VTK file.
        float_format: str
            Float formating.
        """
        coors = self.surf_coors
        cs0, cs1 = coors.shape[0:2]

        nquads = (cs0 - 1) * (cs1 - 1)
        quads = nm.zeros((nquads, 4), dtype=nm.int64)
        kk = 0
        for ii in range(cs0 - 1):
            offs = ii * cs1
            for jj in range(cs1 - 1):
                quads[kk] = nm.array([jj + offs,
                                      jj + offs + cs1,
                                      jj + 1 + offs + cs1,
                                      jj + 1 + offs])
                kk += 1

        f = open(filename, 'w')
        f.write('# vtk DataFile Version 2.0\n')
        f.write('spline surface\nASCII\nDATASET POLYDATA\n')

        ff3 = ' '.join([float_format] * 3) + '\n'
        f.write('POINTS %d float\n' % (cs0 * cs1))
        for ii in range(cs0):
            offs = ii * cs1
            for jj in range(cs1):
                f.write(ff3 % tuple(coors[ii,jj,:]))

        f.write('POLYGONS %d %d\n' % (nquads, nquads * 5))
        for ii in quads:
            f.write('4 %s\n' % (' '.join([str(jj) for jj in ii])))

        f.close()

    def write_control_polygon_vtk(self, filename, float_format='%.6f'):
        """
        Write the control polygon to VTK file.

        Parameters
        ----------
        filename: str
            Name of the VTK file.
        float_format: str
            Float formating.
        """
        coors = self.cp_coors
        cs0, cs1 = coors.shape[0:2]

        lines = []
        nlines = cs0 + cs1
        nlpts = 0
        for ii in range(cs0):
            lines.append(nm.arange(cs1) + ii * cs1)
            nlpts += cs1
        for ii in range(cs1):
            lines.append(nm.arange(cs0) * cs1 + ii)
            nlpts += cs0

        f = open(filename, 'w')
        f.write('# vtk DataFile Version 2.0\n')
        f.write('spline control polygon\nASCII\nDATASET POLYDATA\n')

        ff3 = ' '.join([float_format] * 3) + '\n'
        f.write('POINTS %d float\n' % (cs0 * cs1))
        for ii in range(cs0):
            for jj in range(cs1):
                f.write(ff3 % tuple(coors[ii,jj,:]))

        f.write('LINES %d %d\n' % (nlines, nlines + nlpts))
        for ii in lines:
            f.write('%d %s\n' % (len(ii), ' '.join([str(jj) for jj in ii])))

        f.close()

def get_2d_points(is3d=False):
    """
    Returns the set of points.

    Parameters
    ----------
    is3d : bool
        3D coordinates?
    """
    out = nm.array([(-0.87, 0.15, 0),
                    (-0.70, 0.54, 0),
                    (-0.32, 0.80, 0),
                    (0.15, 0.70, 0),
                    (0.37, 0.26, 0),
                    (0.70, -0.07, 0),
                    (0.67, -0.49, 0),
                    (0.07, -0.81, 0),
                    (-0.44, -0.72, 0),
                    (-0.80, -0.34, 0)])

    if is3d:
        return out

    else:
        return out[:,:2]

def approximation_example():
    """
    The example of using BSplineSurf for approximation
    of the surface given by the set of points.
    """
    # define sample points in 2D
    spts0 = get_2d_points(is3d=True)
    # define sample points in 3D
    spts = nm.array([spts0,
                     spts0 * 0.7 + nm.array([0.1,0,0.5]),
                     spts0 * 0.8 + nm.array([0.2,0,1.5]),
                     spts0 * 1.2 + nm.array([0.4,0,2.5])])

    cyclic=(False, True)
    spl1 = BSplineSurf((3, 3), is_cyclic=cyclic)
    spl1.approximate(spts, (4,8))
    cp = spl1.get_control_points()

    spls2 = BSplineSurf((3, 3), is_cyclic=cyclic)
    spls2.set_control_points(cp)
    spls2.make_knot_vector()
    spls2.set_param_n((12, 24))
    spls2.eval()
    spls2.draw()

def simple_example():
    """
    The example of using B-spline class.
    """
    # define control points in 2D
    cp = get_2d_points()

    spl = BSpline(3, is_cyclic=True)
    spl.set_control_points(cp)
    spl.make_knot_vector()
    spl.set_param_n(150)

    spl.insert_knot(0.7)
    spl.insert_knot(0.7)
    spl.insert_knot(0.7)
    spl.eval()
    spl.draw()

def main(argv):
    # simple_example()
    approximation_example()

if __name__ == '__main__':
    main(sys.argv)

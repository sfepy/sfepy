import numpy as nm
from sfepy.base.base import Struct

class SplineBox(Struct):
    """
    B-spline geometry parametrization. Geometry can be modified
    by moving spline control points.
    """

    @staticmethod
    def mmax(x, y):
        n = len(x)
        aux = nm.zeros((2,n), dtype=nm.int)
        aux[0,:] = x[:]
        aux[1,:] = nm.ones((1,n)) * y

        out = nm.max(aux, axis=0)
        return out

    @staticmethod
    def spsorted(meshsites, sites):
        n1 = len(meshsites)
        n2 = len(sites)
        aux = nm.zeros((n1 + n2, ), dtype=nm.double)
        aux[:n1] = meshsites
        aux[n1:] = sites

        inx = nm.argsort(aux)
        out = nm.where(inx >= len(meshsites)) - nm.arange(len(sites))

        return out[0]

    @staticmethod
    def augknt(knots, k, mults=1):

        if mults > 1:
            aux = []
            for j in k[1:-1]:
                aux += [j] * mults;

        else:
            aux = knots[1:-1]

        augknot = [knots[0]] * k + list(aux) + [knots[-1]] * k

        return augknot

    @staticmethod
    def spcol(knots, k, tau):
        npk = len(knots)
        n = npk - k
        nrows = tau.shape[0]
        pts = tau
        km1 = k - 1

        savl = SplineBox.mmax(SplineBox.spsorted(knots[:n], pts), k)
        b = nm.zeros((nrows, k), dtype=nm.double)
        b[:,0] = nm.ones((nrows,), dtype=nm.double)

        for j in range(km1):
            saved = nm.zeros((nrows,), dtype=nm.double);
            for r in range(j+1):
                tr = knots[savl + r] - pts
                tl = pts - knots[savl + r - j - 1]
                term = nm.double(b[:,r]) / (tr + tl)
                b[:,r] = saved + tr * term
                saved = tl * term
            b[:,j+1] = saved

        idx = nm.where((tau < knots[0]) | (tau > knots[npk-1]))[0]

        if len(idx) > 0:
            b[idx,:] = 0

        idx1 = nm.tile(nm.arange(1-nrows, 1), (k, ))
        idx2 = nm.tile(nrows * savl, (k, ))
        idx3 = nm.tile(nrows * nm.arange(-km1, 1), (nrows, 1))
        idx3 = nm.reshape(idx3, (nrows * k, ), order='F')

        width = n + km1 + km1
        nn = nrows * width
        cc = nm.zeros((nn, ), dtype=nm.double)
        idx = idx1 + idx2 + idx3 - 1
        cc[idx] = b.reshape((len(idx), ), order='F')

        idx1 = nm.tile(nm.arange(1-nrows, 1), (n, ))
        idx2 = nm.tile(nrows * nm.arange(1, n + 1), (nrows, 1))
        idx2 = nm.reshape(idx2, (nrows * n, ), order='F')
        idx = idx1 + idx2 - 1
        colloc = cc[idx].reshape((nrows, n), order='F')

        return colloc

    @staticmethod
    def create_spb(bbox, coors, nsg=None):
        if type(bbox) is not nm.ndarray:
            bbox = nm.array(bbox)

        if type(coors) is not nm.ndarray:
            coors = nm.array(coors)

        dim = coors.shape[1]
        axes =  []

        if nsg is None:
            nsg = nm.ones((dim,), dtype=nm.int)
        else:
            nsg = nm.array(nsg)

        cpt = []
        ncpt = []
        for idim in range(dim):
            axes.append({})
            aux = nm.linspace(bbox[idim,0], bbox[idim,1], nsg[idim] + 1)
            knots = nm.array(SplineBox.augknt(aux, 4))
            axes[idim]['knots'] = knots

            nn = 4 + nsg[idim] - 1
            ncpt.append(nn)
            cpt.append(nm.zeros((nn, ), dtype=nm.double))
            for j in range(nn):
                cpt[idim][j] = nm.sum(knots[(j+1):(j+4)]) / 3.0

            inx = nm.argsort(coors[:,idim])
            aux = SplineBox.spcol(knots, 4, coors[inx,idim])
            axes[idim]['bsc'] = nm.zeros_like(aux)
            axes[idim]['bsc'][inx,:] = aux[:,:]

        ncpts = nm.prod(ncpt)
        cpoints = nm.zeros((ncpts, dim), dtype=nm.double)
        cpoints_idx = nm.zeros(ncpt, dtype=nm.int)
        n2 = nm.prod(ncpt[:2])
        aux = nm.arange(n2).reshape(ncpt[:2], order='F')
        if dim == 2:
            idx = nm.mgrid[0:ncpt[1],0:ncpt[0]]
            cpoints[:,0] = cpt[0][idx[1].reshape((ncpts, ))]
            cpoints[:,1] = cpt[1][idx[0].reshape((ncpts, ))]
            cpoints_idx[:,:] = aux
        elif dim == 3:
            idx = nm.mgrid[0:ncpt[2],0:ncpt[1],0:ncpt[0]]
            cpoints[:,0] = cpt[0][idx[2].reshape((ncpts, ))]
            cpoints[:,1] = cpt[1][idx[1].reshape((ncpts, ))]
            cpoints[:,2] = cpt[2][idx[0].reshape((ncpts, ))]
            for j in range(ncpt[2]):
                cpoints_idx[:,:,j] = aux + j * n2

        return axes, cpoints, cpoints_idx

    def __init__(self, bbox, coors,
                 name='spbox', **kwargs):
        """Create a SplineBox.

        Parameters
        ----------
        bbox : array
            Mesh bounding box.
        coors : array
            Coordinates of mesh nodes.
        name : str
            Object name.
        """
        Struct.__init__(self, name=name, **kwargs)

        self.axes, self.control_points, self.control_points_idx\
                       = self.create_spb(bbox, coors)
        self.dim = self.control_points.shape[1]

        self.control_points0 = self.control_points.copy()

    def get_control_points(self, init=False):
        """Get spline control points coordinates.

        Return
        ------
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
        """Set spline control points position.

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
            self.control_points = cpt_coors.copy()

    def change_shape(self, cpoint, val):
        """Change shape of spline parametrization.

        Parameters
        ----------
        cpoint : list
            The indices of the spline control point.
        val : array
            Displacement.

        """

        idx = self.control_points_idx[cpoint]
        self.control_points[idx] += val

    def evaluate(self, cp_coors=None):
        """Evaluate SplineBox.

        Returns
        -------
        coors : array
            The coordinates corresponding to the actual spline control points position.
        cp_coors : array
            If is not None, use as control points cooardinates.
        """

        coors = nm.zeros((self.axes[0]['bsc'].shape[0], self.dim), dtype=nm.double)

        if cp_coors is None:
            cp_coors = self.control_points

        cptidx = self.control_points_idx
        nn = self.control_points_idx.shape
        if self.dim == 2:
            for i in range(nn[0]):
                for j in range(nn[1]):
                    aux = self.axes[0]['bsc'][:,i] * self.axes[1]['bsc'][:,j]
                    inx = cptidx[i,j]
                    coors += nm.dot(aux[:,nm.newaxis], cp_coors[inx,:][nm.newaxis,:])
        elif self.dim == 3:
            for i in range(nn[0]):
                for j in range(nn[1]):
                    aux = self.axes[0]['bsc'][:,i] * self.axes[1]['bsc'][:,j]
                    for k in range(nn[2]):
                        inx = cptidx[i,j,k]
                        aux2 = aux * self.axes[2]['bsc'][:,k]
                        coors += nm.dot(aux2[:,nm.newaxis], cp_coors[inx,:][nm.newaxis,:])

        return coors

    def dvelocity(self, cpoint, dir):
        """Evaluate derivative of spline in a given control point and direction.

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

        if type(dir) is not nm.ndarray:
            dir = nm.array(dir)

        dvel = nm.zeros((self.axes[0]['bsc'].shape[0], self.dim), dtype=nm.double)

        ax = self.axes
        if self.dim == 2:
            aux = ax[0]['bsc'][:,cpoint[0]] * ax[1]['bsc'][:,cpoint[1]]
            dvel = nm.dot(aux[:,nm.newaxis], dir[nm.newaxis,:])
        elif self.dim == 3:
            aux = ax[0]['bsc'][:,cpoint[0]] * ax[1]['bsc'][:,cpoint[1]]\
                * ax[2]['bsc'][:,cpoint[2]]
            dvel = nm.dot(aux[:,nm.newaxis], dir[nm.newaxis,:])

        return dvel

    def write_vtk(self, filename):

        cptidx = self.control_points_idx
        ncpt = cptidx.shape
        nnd = nm.prod(ncpt)

        f = open(filename, 'w')
        f.write("# vtk DataFile Version 2.6\nspbox file\nASCII\nDATASET UNSTRUCTURED_GRID\n\n")
        f.write("POINTS %d float\n" % nnd)

        if self.dim == 2:

            nel = (ncpt[0] - 1) * ncpt[1] + (ncpt[1] - 1) * ncpt[0]

            for cpt in self.control_points:
                f.write("%e %e 0.0\n" % (cpt[0], cpt[1]))
            f.write("\nCELLS %d %d\n" % (nel, 3 * nel))

            for i in range(ncpt[0]):
                for j in range(ncpt[1] - 1):
                    inx1 = cptidx[i,j]
                    inx2 = cptidx[i,j + 1]
                    f.write("2 %d %d\n" % (inx1, inx2))
            for i in range(ncpt[0] - 1):
                for j in range(ncpt[1]):
                    inx1 = cptidx[i,j]
                    inx2 = cptidx[i + 1,j]
                    f.write("2 %d %d\n" % (inx1, inx2))

        elif self.dim == 3:

            nel = ((ncpt[0] - 1) * ncpt[1] + (ncpt[1] - 1) * ncpt[0]) * ncpt[2]
            nel += ncpt[0] * ncpt[1] * (ncpt[2] - 1)

            for cpt in self.control_points:
                f.write("%e %e %e\n" % (cpt[0], cpt[1], cpt[2]))
            f.write("\nCELLS %d %d\n" % (nel, 3 * nel))

            for k in range(ncpt[2]):
                for i in range(ncpt[0]):
                    for j in range(ncpt[1] - 1):
                        inx1 = cptidx[i, j, k]
                        inx2 = cptidx[i, j + 1, k]
                        f.write("2 %d %d\n" % (inx1, inx2))
                for i in range(ncpt[0] - 1):
                    for j in range(ncpt[1]):
                        inx1 = cptidx[i, j, k]
                        inx2 = cptidx[i + 1, j, k]
                        f.write("2 %d %d\n" % (inx1, inx2))

            for k in range(ncpt[2] - 1):
                for i in range(ncpt[0]):
                    for j in range(ncpt[1]):
                        inx1 = cptidx[i, j, k]
                        inx2 = cptidx[i, j, k + 1]
                        f.write("2 %d %d\n" % (inx1, inx2))

            f.write("\nCELL_TYPES %d\n" % nel )
            f.write("3\n" * nel)
            f.close

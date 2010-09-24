import numpy as nm
import scipy.sparse as sp

from sfepy.base.base import Struct, assert_
from sfepy.base.la import permutations
import extmods.meshutils as mu

def _build_orientation_map(n_fp):
    """
    The keys are binary masks of the lexicographical ordering of facet
    vertices. A bit i set to one means `v[i] < v[i+1]`.

    The values are `[original_order, permutation, orientation]`, where
    `permutation` can be used to sort facet vertices lexicographically,
    and `orientation` is the order of the first vertex + 1 times the
    sign. Hence `permuted_facet = facet[permutation]`.
    """
    indices = range(n_fp)

    cmps = [(i1, i2) for i2 in indices for i1 in indices[:i2]]
    powers = [2**ii for ii in range(len(cmps))]

    ori_map = {}
    for indx in permutations(indices):
        key = 0
        sign = 1
        for ip, power in enumerate(powers):
            i1, i2 = cmps[ip]
            less = (indx[i1] < indx[i2])
            key += power * less
            if not less:
                sign *= -1

        isort = nm.argsort(indx)
        ori_map[key] = [indx, isort, sign * (indx[0] + 1)]

    return ori_map, cmps, powers

def _get_signed_orientation(ori, ori_map):
    """
    Transform orientation according to `ori_map`, i.e. from bit mask
    encoding to signed values.
    """
    i_from = nm.array(ori_map.keys())
    i_to = [ii[-1] for ii in ori_map.values()]

    i_map = nm.zeros((i_from.max() + 1,), dtype=nm.int8)
    i_map[i_from] = i_to

    signed_ori = i_map[ori]
    assert_((signed_ori != 0).all())

    return signed_ori

def _permute_facets(facets, ori, ori_map):
    """
    Return a copy of `facets` array with vertices sorted lexicographically.
    """
    assert_((nm.setmember1d(nm.unique(ori), ori_map.keys())).all())

    permuted_facets = facets.copy()

    for key, ori_map in ori_map.iteritems():
        perm = ori_map[1]
        ip = nm.where(ori == key)[0]
        for ic0, ic1 in enumerate(perm):
            permuted_facets[ip,ic0] = facets[ip,ic1]

    return permuted_facets

def _orient_facets(ori, facets, cmps, powers):
    for ip, power in enumerate(powers):
        i1, i2 = cmps[ip]
        iw = nm.where(facets[:,i1] < facets[:,i2])[0]
        ori[iw] += power

class Facets(Struct):

    @staticmethod
    def from_domain(domain, kind):
        groups = domain.groups

        if kind == 'edges':
            n_obj = [group.shape.n_edge_total for group in groups.itervalues()]
            nn = 2
        else:
            n_obj = [group.shape.n_face_total for group in groups.itervalues()]
            nn = 4

        n_all_obj = sum(n_obj)
        indices = nm.zeros((n_all_obj, 3), nm.int32)
        all_facets = nm.empty((n_all_obj, nn), nm.int32)
        all_facets.fill(-1)

        single_facets = {}

        ii = 0
        for ig, group in groups.iteritems():
            conn, gel = group.conn, group.gel
            n_el = group.shape.n_el
            n_all_item = n_obj[ig]

            if kind == 'edges':
                n_item, items = gel.n_edge, gel.edges

            else:
                n_item, items = gel.n_face, gel.faces

            n_fp = items.shape[1]
            single_facets[ig] = items

            io = slice(ii, ii + n_all_item)

            indices[io,0] = ig

            ie = nm.arange(n_el, dtype=nm.int32)
            ie = nm.repeat(ie, n_item)
            indices[io,1] = ie

            iobj = nm.arange(n_item, dtype=nm.int32)
            indices[io,2] = nm.tile(iobj, n_el)

            facets = conn[:, items]
            facets = facets.reshape((n_all_item, n_fp))

            all_facets[io,:n_fp] = facets

            ii += n_all_item

        if (ii != sum( n_obj )):
            msg = 'neighbour_list size mismatch! (%d == %d = sum( %s ))'\
                  % (ii, sum( n_obj ), n_obj)
            raise ValueError( msg )

        obj = Facets('facets', kind, domain, single_facets,
                     n_obj, indices, all_facets)

        return obj

    def __init__(self, name, kind, domain, single_facets,
                 n_obj, indices, facets):
        Struct.__init__(self, name=name, kind=kind, domain=domain,
                        single_facets=single_facets,
                        n_obj=n_obj, indices=indices, facets=facets)
        self.n_all_obj, self.n_col = facets.shape
        self.n_gr = len(self.n_obj)

        self.indx = {}
        ii = 0
        for ig, nn in enumerate(self.n_obj):
            self.indx[ig] = slice(ii, ii+nn)
            ii += nn

        self.n_fps_vec = nm.empty(self.n_all_obj, dtype=nm.int32)
        self.n_fps = {}
        for ig, facet in self.single_facets.iteritems():
            self.n_fps_vec[self.indx[ig]] = facet.shape[1]
            self.n_fps[ig] = facet.shape[1]

    def sort_and_orient(self):
        all_permuted_facets = nm.empty((self.n_all_obj + 2, self.n_col),
                                     dtype=nm.int32)
        all_permuted_facets.fill(-1)

        sentinel = self.domain.shape.n_nod

        aux = nm.repeat(nm.array([sentinel], nm.int32), self.n_col)
        all_permuted_facets[-2] = aux
        all_permuted_facets[-1] = aux + 1
        oris = {}
        ori_maps = {}

        for ig in range(self.n_gr):
            io = self.indx[ig]
            facets = self.facets[io]

            n_fp = self.n_fps[ig]
            ori_map, cmps, powers = _build_orientation_map(n_fp)

            # Determine orientation.
            ori = nm.zeros((facets.shape[0],), dtype=nm.int8)
            _orient_facets(ori, facets, cmps, powers)

            # Permute each facet to have indices in ascending order, so
            # that lexicographic sorting works.
            permuted_facets = _permute_facets(facets, ori, ori_map)
            all_permuted_facets[io] = permuted_facets

            ori = _get_signed_orientation(ori, ori_map)

            oris[ig] = ori
            ori_maps[ig] = ori_map

        self.permuted_facets = all_permuted_facets
        self.oris = oris
        self.ori_maps = ori_maps

    def setup_unique(self):
        """
        `sorted_facets` == `permuted_facets[perm]`
        `permuted_facets` == `sorted_facets[perm_i]`
        `uid` : unique id in order of `sorted_facets`
        `uid_i` : unique id in order of `permuted_facets` or `facets`
        """
        ii = nm.arange(self.permuted_facets.shape[0], dtype=nm.int32)
        aux = nm.concatenate((self.permuted_facets, ii[:,None]), 1).copy()
        sort_cols = nm.arange(self.permuted_facets.shape[1],
                              dtype=nm.int32)

        mu.sort_rows(aux, sort_cols)

        self.perm = perm = aux[:,-1].copy()
        self.sorted_facets = aux[:,:-1].copy()
        aux = nm.arange(perm.shape[0], dtype=nm.int32)
        self.perm_i = nm.zeros_like(self.perm)
        self.perm_i[perm] = aux

        ic = nm.where(nm.abs(nm.diff(self.sorted_facets, axis=0)).sum(1), 0, 1)
        ic = ic.astype(nm.int32)

        self.n_unique = len(ic) - nm.sum(ic) - 1
        self.unique_list = nm.where(ic[:-1] == 0)[0].astype(nm.int32)
        assert_(len(self.unique_list) == self.n_unique)

        ii = nm.cumsum( ic[:-1] == 0, dtype = nm.int32 )
        self.uid = ii.copy()
        self.uid[0], self.uid[1:] = 0, ii[:-1]
        self.uid_i = self.uid[self.perm_i[:-2]]

    def setup_neighbours(self):
        """
        For each unique facet:
           - indices of facets - sparse matrix (n_unique x n_all_obj)
             mtx[i, j] == 1 if facet[j] has uid[i]
           - number of elements it is in
        """
        ones = nm.ones((self.n_all_obj,), dtype=nm.bool)

        self.mtx = sp.coo_matrix((ones, (self.uid, self.perm[:-2])))

        self.n_in_el = self.mtx * ones.astype(nm.int32)

    def mark_surface_facets(self):
        """
        flag: 0 .. inner, 2 .. edge, 3 .. triangle, 4 .. quadrangle
        """
        nn = self.n_in_el[self.uid_i]

        flag = self.n_fps_vec * (nn == 1)

        return flag

    def get_coors(self, ig=0):
        """
        Get the coordinates of vertices of unique facets in group `ig`.

        Returns
        -------
        coors : array
            The coordinates in an array of shape `(n_f, n_v, dim)`.
        uid : array
            The unique ids of facets in the order of `coors`.
        """
        cc = self.domain.get_mesh_coors()

        uid_i = self.uid_i[self.indx[ig]]
        uid, ii = nm.unique1d(uid_i, return_index=True)

        coors = cc[self.facets[ii,:self.n_fps[ig]]]

        return coors, uid

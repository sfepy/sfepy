import numpy as nm
import scipy.sparse as sp

from sfepy.base.base import Struct, dict_to_array, assert_
from sfepy.base.compat import in1d, unique
from sfepy.linalg import permutations, map_permutations, insert_strided_axis

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
    assert_((in1d(nm.unique(ori), ori_map.keys())).all())

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

_quad_ori_groups = {
    0 : 0,
    1 : 0,
    3 : 7,
    4 : 0,
    6 : 7,
    7 : 7,
    11 : 11,
    15 : 11,
    20 : 52,
    22 : 30,
    30 : 30,
    31 : 63,
    32 : 33,
    33 : 33,
    41 : 33,
    43 : 11,
    48 : 56,
    52 : 52,
    56 : 56,
    57 : 56,
    59 : 63,
    60 : 52,
    62 : 30,
    63 : 63,
}

_quad_orientations = {
    0  : [0, 1, 3, 2],
    7  : [3, 2, 0, 1],
    11 : [2, 3, 0, 1],
    30 : [1, 0, 2, 3],
    33 : [1, 0, 3, 2],
    52 : [2, 3, 1, 0],
    56 : [3, 2, 1, 0],
    63 : [0, 1, 2, 3],
}

def get_facet_dof_permutations(int_coors, ori_maps):
    """
    Prepare DOF permutation vector for each possible facet orientation.
    """
    dof_perms = {}

    ori = nm.empty((int_coors.shape[0],), dtype=nm.int32)
    for ig, ori_map in ori_maps.iteritems():
        dof_perms[ig] = {}
        for key in ori_map.iterkeys():
            ori.fill(key)
            permuted_int_coors = _permute_facets(int_coors, ori, ori_map)

            perm = map_permutations(permuted_int_coors, int_coors,
                                    check_same_items=True)

            dof_perms[ig][key] = perm

    return dof_perms

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
            if n_obj[ig] == 0:
                single_facets[ig] = nm.array([[]], dtype=nm.int32)
                continue

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
        signed_oris = {}
        ori_maps = {}

        for ig in range(self.n_gr):
            if self.n_obj[ig] == 0: continue
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

            signed_ori = _get_signed_orientation(ori, ori_map)

            oris[ig] = ori
            signed_oris[ig] = signed_ori
            ori_maps[ig] = ori_map

        self.permuted_facets = all_permuted_facets
        self.oris = oris
        self.signed_oris = signed_oris
        self.ori_maps = ori_maps

    def setup_unique(self):
        """
        `sorted_facets` == `permuted_facets[perm]`
        `permuted_facets` == `sorted_facets[perm_i]`
        `uid` : unique id in order of `sorted_facets`
        `uid_i` : unique id in order of `permuted_facets` or `facets`
        """
        self.perm = nm.lexsort(self.permuted_facets.T[::-1])
        self.sorted_facets = self.permuted_facets[self.perm]

        self.perm_i = nm.zeros_like(self.perm)
        self.perm_i[self.perm] = nm.arange(self.perm.shape[0], dtype=nm.int32)

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

    def find_group_interfaces(self, return_surface=True):
        """
        Find facets that create boundary between different element
        groups, i.e. facets that each belongs to two elements in
        different groups.

        Parameters
        ----------
        return_surface : bool
            If True, the surface facets are also returned.

        Returns
        -------
        inter_facets : array
            The array with indices to `self.facets` of shape `(n_i, 2)`,
            where `n_i` is the number of the interface facets. Each row
            corresponds to a single unique facet, each column to the
            corresponding two facets from each side of the interface.
        surface_facets : array, optional
            The array with indices to `self.facets` of shape `(n_s,)`,
            where `n_s` is the number of the surface facets.
        """
        mtx = self.mtx.tocsr()

        # ... inner facets are in two elements
        i2 = nm.where(self.n_in_el == 2)[0]

        face_map = mtx[i2]
        uid, ifacets = face_map.nonzero()

        # ... sort by uid
        ii = nm.argsort(uid)
        ifacets = ifacets[ii]
        ifacets.shape = (i2.shape[0], 2)

        igs = self.indices[ifacets, 0]

        # ... interface facets are in two groups
        ii = nm.where(igs[:, 0] != igs[:, 1])[0]
        inter_facets = ifacets[ii]

        out = [inter_facets]

        if return_surface:
            i1 = nm.where(self.n_in_el == 1)[0]
            face_map = mtx[i1]
            _, surface_facets = face_map.nonzero()
            out = out + [surface_facets]

        return out

    def mark_surface_facets(self):
        """
        flag: 0 .. inner, 2 .. edge, 3 .. triangle, 4 .. quadrangle
        """
        nn = self.n_in_el[self.uid_i]

        flag = self.n_fps_vec * (nn == 1)

        return flag

    def get_coors(self, ig=None):
        """
        Get the coordinates of vertices of unique facets in group `ig`.

        Parameters
        ----------
        ig : int, optional
            The element group. If None, the coordinates for all groups
            are returned, filled with zeros at places of missing
            vertices, i.e. where facets having less then the full number
            of vertices (`n_v`) are.

        Returns
        -------
        coors : array
            The coordinates in an array of shape `(n_f, n_v, dim)`.
        uid : array
            The unique ids of facets in the order of `coors`.
        """
        cc = self.domain.get_mesh_coors()

        if ig is None:
            uid, ii = unique(self.uid_i, return_index=True)
            facets = self.facets[ii]
            aux = insert_strided_axis(facets, 2, cc.shape[1])
            coors = nm.where(aux >= 0, cc[facets], 0.0)

        else:
            uid_i = self.uid_i[self.indx[ig]]
            uid, ii = unique(uid_i, return_index=True)

            coors = cc[self.facets[ii,:self.n_fps[ig]]]

        return coors, uid

    def get_uid_per_elements(self, ig=0):
        """
        Get unique ids of facets for all elements in group `ig`.
        """
        n_facet = self.single_facets[ig].shape[0]

        uid_i = self.uid_i[self.indx[ig]]
        uid_i.shape = (uid_i.shape[0] / n_facet, n_facet)

        return uid_i

    def get_dof_orientation_maps(self, nodes):
        """
        Given description of facet DOF nodes, return the corresponding
        integer coordinates and orientation maps.

        Notes
        -----
        Assumes single facet type in all groups.
        """
        inod = nm.arange(self.n_fps[0], dtype=nm.int32)

        int_coors = nodes[0][:, inod]

        if int_coors.shape[1] <= 3: # Simplex facet.
            ori_maps = self.ori_maps

        else: # Tensor product facet.
            ori_maps = {}
            for ig, ori_map in self.ori_maps.iteritems():
                ori_maps[ig] = {}
                for key in ori_map.iterkeys():
                    new_key = _quad_ori_groups[key]
                    ori_maps[ig][key] = [None, _quad_orientations[new_key]]

        return int_coors, ori_maps

    def get_facet_dof_permutations(self, nodes):
        """
        Given description of facet DOF nodes, return the DOF
        permutations for all possible facet orientations.
        """
        int_coors, ori_maps = self.get_dof_orientation_maps(nodes)

        aux = get_facet_dof_permutations(int_coors, ori_maps)
        dof_perms = {}
        for ig, val in aux.iteritems():
            dof_perms[ig] = dict_to_array(val)

        return dof_perms

    def get_complete_facets(self, vertices, ig=0, mask=None):
        """
        Get complete facets in group `ig` that are defined by the given
        vertices, or mask, if given.

        Parameters
        ----------
        vertices : array
            The list of vertices.
        ig : int
            The group index.
        mask : array, optional
            Alternatively to `vertices`, a mask can be given with 1 at
            indices equal to the vertices and 0 elsewhere.

        Returns
        -------
        ifacets : array
            The indices into `self.facets`.
        """
        n_fp = self.n_fps[ig]
        indx = self.indx[ig]

        ii = self.facets[indx,:n_fp]

        if mask is None:
            mask = nm.zeros(ii.max()+1, dtype=nm.bool)
            mask[vertices] = True

        aux = nm.sum(mask[ii], 1)

        # Points to self.facets.
        ifacets = indx.start + nm.where(aux == n_fp)[0]

        return ifacets

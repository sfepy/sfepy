from sfepy.base.base import *

class FESurface(Struct):
    """Description of a surface of a finite element domain."""

    def __init__(self, name, region, efaces, volume_econn, ig):
        """nodes[leconn] == econn"""
        """nodes are sorted by node number -> same order as region.vertices"""
        self.name = get_default(name, 'surface_data_%s' % region.name)

        face_indices = region.fis[ig]
        
        faces = efaces[face_indices[:,1]]
        if faces.size == 0:
            raise ValueError('region with group with no faces! (%s)'
                             % region.name)

        try:
            ee = volume_econn[face_indices[:,0]]
        except:
            pdb.set_trace()

        econn = nm.empty(faces.shape, dtype=nm.int32)
        for ir, face in enumerate( faces ):
#            print ir, face, ee[ir,face]
            econn[ir] = ee[ir,face]

        ef = econn.flat
        nodes = nm.unique1d(ef)

        aux = -nm.ones((nm.max( ef ) + 1,), dtype=nm.int32)
        aux[nodes] = nm.arange(len(nodes), dtype=nm.int32)
        leconn = aux[econn].copy()
        assert_(nm.alltrue(nodes[leconn] == econn))
        
        n_fa, n_fp = face_indices.shape[0], faces.shape[1]
        face_type = 's%d' % n_fp

        # Store bkey in SurfaceData, so that base function can be
        # queried later.
        bkey = 'b%s' % face_type[1:]

        self.econn = econn
        self.fis = face_indices
        self.n_fa, self.n_fp = n_fa, n_fp
        self.nodes = nodes
        self.leconn = leconn
        self.face_type = face_type
        self.bkey = bkey

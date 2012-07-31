"""
GeometryElement describes the geometric entities of a finite element mesh.

Notes
-----
* geometry_data: surface facets are assumed to be of the same kind for
  each geometry element - wedges or pyramides are not supported.
* the orientation is a tuple:
  (root1, vertices of direction vectors, swap from, swap to, root2, ...)
"""
import numpy as nm

from sfepy.base.base import assert_, Struct

geometry_data = {
    '1_2' : Struct(coors = [[0.0],
                            [1.0]],
                   conn = [0, 1],
                   faces = None,
                   edges = None,
                   volume = 1.0,
                   orientation = None,
                   surface_facet_name = None),

    '2_3' : Struct(coors = [[0.0, 0.0],
                            [1.0, 0.0],
                            [0.0, 1.0]],
                   conn = [0, 1, 2],
                   faces = None,
                   edges = [[0, 1],
                            [1, 2],
                            [2, 0]],
                   volume = 0.5,
                   orientation = (0, (1, 2), 1, 2),
                   surface_facet_name = '1_2'),

    '2_4' : Struct(coors = [[0.0, 0.0],
                            [1.0, 0.0],
                            [1.0, 1.0],
                            [0.0, 1.0]],
                   conn = [0, 1, 2, 3],
                   faces = None,
                   edges = [[0, 1],
                            [1, 2],
                            [2, 3],
                            [3, 0]],
                   volume = 1.0,
                   # Not finished...
                   orientation = (0, (1, 3), (0, 1), (3, 2)),
                   surface_facet_name = '1_2'),

    '3_4' : Struct(coors = [[0.0, 0.0, 0.0],
                            [1.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0],
                            [0.0, 0.0, 1.0]],
                   conn = [0, 1, 2, 3],
                   faces = [[0, 2, 1],
                            [0, 3, 2],
                            [0, 1, 3],
                            [1, 2, 3]],
                   edges = [[0, 1],
                            [1, 2],
                            [2, 0],
                            [0, 3],
                            [1, 3],
                            [2, 3]],
                   volume = 1.0 / 6.0,
                   orientation = (0, (1, 2, 3), 0, 3),
                   surface_facet_name = '2_3'),

    '3_8' : Struct(coors = [[0.0, 0.0, 0.0],
                            [1.0, 0.0, 0.0],
                            [1.0, 1.0, 0.0],
                            [0.0, 1.0, 0.0],
                            [0.0, 0.0, 1.0],
                            [1.0, 0.0, 1.0],
                            [1.0, 1.0, 1.0],
                            [0.0, 1.0, 1.0]],
                   conn = [0, 1, 2, 3, 4, 5, 6, 7],
                   faces = [[0, 3, 2, 1],
                            [0, 4, 7, 3],
                            [0, 1, 5, 4],
                            [4, 5, 6, 7],
                            [1, 2, 6, 5],
                            [3, 7, 6, 2]],
                   edges = [[0, 1],
                            [1, 2],
                            [2, 3],
                            [3, 0],
                            [4, 5],
                            [5, 6],
                            [6, 7],
                            [7, 4],
                            [0, 4],
                            [1, 5],
                            [2, 6],
                            [3, 7]],
                   volume = 1.0,
                   # Not finished...
                   orientation = (0, (1, 3, 4), (0, 1, 2, 3), (4, 5, 6, 7) ),
                   surface_facet_name = '2_4'),
}

def setup_orientation(vecs_tuple):
    cycle = range(len(vecs_tuple) / 4)

    roots = nm.array([vecs_tuple[4*ii] for ii in cycle], dtype=nm.int32)
    vecs = nm.array([vecs_tuple[4*ii+1] for ii in cycle],
                    dtype=nm.int32, ndmin=2)
    swap_from = nm.array([vecs_tuple[4*ii+2] for ii in cycle],
                         dtype=nm.int32, ndmin=2)
    swap_to = nm.array([vecs_tuple[4*ii+3] for ii in cycle],
                       dtype=nm.int32, ndmin=2)

    return roots, vecs, swap_from, swap_to

class GeometryElement(Struct):
    """
    The geometric entities of a finite element mesh.
    """

    def __init__(self, name):
        """
        Parameters
        ----------
        name : str
            The name of the entity, one of the keys in geometry_data
            dictionary.
        """
        self.name = name

        gd = geometry_data[name]

        self.coors = nm.array(gd.coors, dtype=nm.float64)

        self.conn = nm.array(gd.conn, dtype=nm.int32)

        self.n_vertex, self.dim = self.coors.shape
        self.is_simplex = self.n_vertex == (self.dim + 1)

        self.vertices = nm.arange(self.n_vertex, dtype=nm.int32)

        if gd.edges is not None:
            self.edges = nm.array(gd.edges, dtype=nm.int32)
            self.n_edge = self.edges.shape[0]
        else:
            self.edges = gd.edges
            self.n_edge = 0

        if gd.faces is not None:
            self.faces = nm.array(gd.faces, dtype=nm.int32)
            self.n_face = self.faces.shape[0]
        else:
            self.faces = gd.faces
            self.n_face = 0

        if gd.orientation is not None:
            aux = setup_orientation(gd.orientation)
            self.orientation = Struct(name='orientation',
                                      roots=aux[0], vecs=aux[1],
                                      swap_from=aux[2], swap_to=aux[3])
        else:
            self.orientation = None

        self.surface_facet_name = gd.surface_facet_name
        self.surface_facet = None

    def get_interpolation_name(self):
        """
        Get the name of corresponding linear interpolant.
        """
        if self.is_simplex:
            suffix = '_P1'
        else:
            suffix = '_Q1'
        return self.name + suffix

    def get_surface_entities(self):
        """
        Return self.vertices in 1D, self.edges in 2D and self.faces in 3D.
        """
        if self.dim == 1:
            return self.vertices
        elif self.dim == 2:
            return self.edges
        else:
            assert_(self.dim == 3)
            return self.faces

    def get_edges_per_face(self):
        """
        Return the indices into self.edges per face.
        """
        if self.dim == 3:
            # Assign edges to a face (in order).
            indx = {3: [[0, 1], [1, 2], [2, 0]],
                    4: [[0, 1], [1, 2], [2, 3], [3, 0]]}
            epf = []
            se = [set(edge) for edge in self.edges]
            iis = indx[self.surface_facet.n_vertex]
            for face in self.faces:
                aux = []

                for ii in iis:
                    edge = set(face[ii])
                    ie = se.index(edge)
                    aux.append(ie)

                epf.append(aux)

        else:
            epf = nm.arange(self.edges.shape[0])[:,nm.newaxis]

        return nm.array(epf, dtype=nm.int32)

    def get_conn_permutations(self):
        """
        Get all possible connectivity permutations corresponding to different
        spatial orientations of the geometry element.
        """
        if self.dim < 3:
            perms = [nm.roll(self.conn, -ii) for ii in range(self.n_vertex)]
            perms = nm.vstack(perms)

        else:
            _perms3d = {
                '3_4' : [[0, 1, 2, 3],
                         [1, 2, 0, 3],
                         [2, 0, 1, 3],
                         [1, 3, 2, 0],
                         [3, 0, 2, 1],
                         [3, 1, 0, 2],
                         [2, 1, 3, 0],
                         [0, 2, 3, 1],
                         [0, 3, 1, 2]],
                '3_8' : [[0, 1, 2, 3, 4, 5, 6, 7],
                         [1, 2, 3, 0, 5, 6, 7, 4],
                         [2, 3, 0, 1, 6, 7, 4, 5],
                         [3, 0, 1, 2, 7, 4, 5, 6],
                         [3, 2, 6, 7, 0, 1, 5, 4],
                         [7, 6, 5, 4, 3, 2, 1, 0],
                         [4, 5, 1, 0, 7, 6, 2, 3],
                         [1, 5, 6, 2, 0, 4, 7, 3],
                         [5, 4, 7, 6, 1, 0, 3, 2],
                         [4, 0, 3, 7, 5, 1, 2, 6]],
            }

            perms = nm.array(_perms3d[self.name], dtype=nm.int32)

        return perms

    def create_surface_facet(self):
        """
        Create a GeometryElement instance corresponding to this instance
        surface facet.
        """
        self.surface_facet = GeometryElement(self.surface_facet_name)

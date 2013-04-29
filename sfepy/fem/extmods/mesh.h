#ifndef _MESH_H_
#define _MESH_H_ 1

#ifdef __cplusplus
#  define BEGIN_C_DECLS         extern "C" {
#  define END_C_DECLS           }
#else
#  define BEGIN_C_DECLS
#  define END_C_DECLS
#endif

#include "types.h"
#include "version.h"

// Inspired by A. Logg: Efficient Representation of Computational Meshes.

// Entity, dimension, codimension. Max. dimension is D.
#define Vertex 0 // 0,     D
#define Edge   1 // 1,     D - 1
#define Face   2 // 2,     D - 2
#define Cell   3 // D,     0
#define Facet  4 // D - 1, 1

typedef struct MeshGeometry {
  uint32 num;
  uint32 dim;
  float64 *coors; // Shape: (num, dim) by rows.
} MeshGeometry;

typedef struct MeshConnectivity {
  uint32 num; // MeshTopology->num[i] - number of entities.
  uint32 n_incident; // Total number of incident entities.
  uint32 *indices; // Length: num.
  uint32 *offsets; // Length: MeshTopology->num[i] + 1.
  uint32 offset; // For homogeneous connectivity (e.g. edges). (?)
} MeshConnectivity;

typedef struct ConnIter {
  uint32 ii;
  uint32 num;
  uint32 *ptr;
  MeshConnectivity *conn;
} ConnIter;

// Connectivity defines incidence.
// conn[IJ[d1, d2]] are entities of dimension d2 incident to entities of
// dimension d1.

// d1 > d2:
// (d2, i2) is incident to (d1, i1) if (d2, i2) is contained in (d1, i1)

// d1 < d2:
// (d2, i2) is incident to (d1, i1) if (d1, i1) is incident to (d2, i2)

// d1 = d2, d1, d2 > 0:
// (d2, i2) is incident to (d1, i1) if both are incident to a common vertex

// d1 = d2 = 0:
// (d2, i2) is incident to (d1, i1) if both are incident to a common cell

#define IJ(D, d1, d2) (D) * (d1) + (d2)
typedef struct MeshTopology {
  uint32 max_dim;
  uint32 num[4]; // Number of entities of given dimension.
  MeshConnectivity _conn[16];
  MeshConnectivity *conn[16];
} MeshTopology;

typedef struct EntityIter {
  uint32 ii;
  uint32 dim;
  MeshTopology *topology;
} EntityIter;

typedef struct Mesh {
  MeshGeometry geometry[1];
  MeshTopology topology[1];
} Mesh;

typedef struct Region {
  uint32 *vertices;
  uint32 *edges;
  uint32 *faces;
  uint32 *cells;
} Region;

int32 mesh_create(Mesh **p_mesh);
int32 mesh_destroy(Mesh **p_mesh);
int32 mesh_init(Mesh *mesh);
int32 mesh_print(Mesh *mesh, FILE *file, int32 header_only);

int32 conn_iter_init(ConnIter *iter, MeshConnectivity *conn);
int32 conn_iter_print(ConnIter *iter, FILE *file);
int32 conn_iter_next(ConnIter *iter);
int32 conn_print(MeshConnectivity *conn, FILE *file);

int32 entity_iter_init(EntityIter *iter, int32 dim, MeshTopology *topology);
int32 entity_iter_print(EntityIter *iter, FILE *file);
int32 entity_iter_next(EntityIter *iter);

int32 mesh_set_coors(Mesh *mesh, float64 *coors, int32 num, int32 dim);

int32 mesh_setup_connectivity(Mesh *mesh, int32 d1, int32 d2);
int32 mesh_build(Mesh *mesh, int32 dim);
int32 mesh_transpose(Mesh *mesh, int32 d1, int32 d2);
int32 mesh_intersect(Mesh *mesh, int32 d1, int32 d2, int32 d3);

// Mesh entity is given by (dimension, index) or (dim, ii).
// Index ii is in [0, num[dim] - 1].
int32 mesh_get_entity(int32 *n_items, int32 dim, int32 ii);

// ?
int32 mesh_iterate_entities(Mesh *mesh, int32 *entity, int32 dim, Region *reg);

#endif /* !MESH_H */

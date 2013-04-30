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

// Special uint32 values meaning "not set".
#define UINT32_None -1

typedef struct Mesh Mesh;

typedef struct MeshEntity {
  uint32 dim; // Topological dimension.
  uint32 ii;  // Local index within entities of the given topological
              // dimension.
  Mesh *mesh;
} MeshEntity;

typedef struct MeshEntityIterator {
  uint32 it; // Current iteration position.
  uint32 it_end; // End iteration position.
  uint32 *ptr; // If given, entity->ii = ptr[it].
  MeshEntity entity[1];
} MeshEntityIterator;

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
} MeshConnectivity;

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

#define IJ(D, d1, d2) ((D + 1)) * (d1) + (d2)
typedef struct MeshTopology {
  uint32 max_dim;
  uint32 num[4]; // Number of entities of given dimension.
  MeshConnectivity _conn[16];
  MeshConnectivity *conn[16];
} MeshTopology;

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

int32 mei_init(MeshEntityIterator *iter, Mesh *mesh, uint32 dim);
int32 mei_init_conn(MeshEntityIterator *iter, MeshEntity *entity, uint32 dim);
int32 mei_print(MeshEntityIterator *iter, FILE *file);
int32 mei_go(MeshEntityIterator *iter);
int32 mei_next(MeshEntityIterator *iter);

int32 conn_alloc(MeshConnectivity *conn, uint32 num, uint32 n_incident);
int32 conn_free(MeshConnectivity *conn);
int32 conn_print(MeshConnectivity *conn, FILE *file);
int32 conn_set_to_free(MeshConnectivity *conn, uint32 ii, uint32 incident);

int32 mesh_set_coors(Mesh *mesh, float64 *coors, int32 num, int32 dim);

int32 mesh_setup_connectivity(Mesh *mesh, int32 d1, int32 d2);
int32 mesh_free_connectivity(Mesh *mesh, int32 d1, int32 d2);
int32 mesh_build(Mesh *mesh, int32 dim);
int32 mesh_transpose(Mesh *mesh, int32 d1, int32 d2);
int32 mesh_intersect(Mesh *mesh, int32 d1, int32 d2, int32 d3);

// Mesh entity is given by (dimension, index) or (dim, ii).
// Index ii is in [0, num[dim] - 1].
int32 mesh_get_entity(int32 *n_items, int32 dim, int32 ii);

// ?
int32 mesh_iterate_entities(Mesh *mesh, int32 *entity, int32 dim, Region *reg);

#endif /* !MESH_H */

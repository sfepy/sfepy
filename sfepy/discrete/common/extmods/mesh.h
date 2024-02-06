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
#define UINT32_None (uint32)-1

// Pointer to indices + number of items.
typedef struct Indices {
  uint32 *indices;
  uint32 num;
} Indices;

// Pointer to mask + number of items.
typedef struct Mask {
  char *mask;
  uint32 num;
  uint32 n_true;
} Mask;

typedef struct MeshGeometry {
  uint32 num;
  uint32 dim;
  float64 *coors; // Shape: (num, dim) by rows.
} MeshGeometry;

typedef struct MeshConnectivity {
  uint32 num; // MeshTopology->num[i] - number of entities.
  uint32 n_incident; // Total number of incident entities.
  uint32 *indices; // Length: n_incident.
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
  uint32 *cell_types; // Cell types = indices to LocalEntities.
  uint32 *face_oris; // Face orientations.
  uint32 *edge_oris; // Edge orientations.
  MeshConnectivity _conn[16];
  MeshConnectivity *conn[16];
} MeshTopology;

#define MAX_EL_TYPES 6
#define Bar 0
#define Triangle 1
#define Quadrilateral 2
#define Tetrahedron 3
#define Hexahedron 4
#define Wedge 5

// Facets for various reference element types.
typedef struct LocalEntities {
  uint32 num; // Lengths of edges and faces = MAX_EL_TYPES.
  MeshConnectivity _edges[MAX_EL_TYPES];
  MeshConnectivity *edges[MAX_EL_TYPES];
  MeshConnectivity _faces[MAX_EL_TYPES];
  MeshConnectivity *faces[MAX_EL_TYPES];
} LocalEntities;

typedef struct Mesh {
  MeshGeometry geometry[1];
  MeshTopology topology[1];
  LocalEntities entities[1];
} Mesh;

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

typedef struct Region {
  Indices *vertices;
  Indices *edges;
  Indices *faces;
  Indices *cells;
} Region;

int32 mesh_create(Mesh **p_mesh);
int32 mesh_destroy(Mesh **p_mesh);
int32 mesh_init(Mesh *mesh);
int32 mesh_free(Mesh *mesh);
int32 mesh_print(Mesh *mesh, FILE *file, int32 header_only);

int32 mei_init(MeshEntityIterator *iter, Mesh *mesh, uint32 dim);
int32 mei_init_conn(MeshEntityIterator *iter, MeshEntity *entity, uint32 dim);
int32 mei_init_sub(MeshEntityIterator *iter, Mesh *mesh,
                   Indices *entities, uint32 dim);
int32 mei_print(MeshEntityIterator *iter, FILE *file);
int32 mei_go(MeshEntityIterator *iter);
int32 mei_next(MeshEntityIterator *iter);

int32 ind_print(Indices *ind, FILE *file);

int32 conn_alloc(MeshConnectivity *conn, uint32 num, uint32 n_incident);
int32 conn_resize(MeshConnectivity *conn, uint32 num, uint32 n_incident);
int32 conn_free(MeshConnectivity *conn);
int32 conn_print(MeshConnectivity *conn, FILE *file);
int32 conn_set_from(MeshConnectivity *conn, MeshConnectivity *other);
int32 conn_set_to_free(MeshConnectivity *conn, uint32 ii, uint32 incident);

int32 mesh_set_coors(Mesh *mesh, float64 *coors, int32 num, int32 dim,
                     int32 tdim);

int32 mesh_setup_connectivity(Mesh *mesh, int32 d1, int32 d2);
int32 mesh_free_connectivity(Mesh *mesh, int32 d1, int32 d2);
int32 mesh_build(Mesh *mesh, int32 dim);
int32 mesh_transpose(Mesh *mesh, int32 d1, int32 d2);
int32 mesh_intersect(Mesh *mesh, int32 d1, int32 d2, int32 d3);

// Count non-unique entities of dimension `dim` incident to `entities` of
// dimension `dent`.
uint32 mesh_count_incident(Mesh *mesh, int32 dim,
                           Indices *entities, int32 dent);

// Get non-unique entities `incident` of dimension `dim` incident to `entities`
// of dimension `dent`. As each of entities can be in several entities of
// dimension `dent`, `incident` is stored in MeshConnectivity structure.
// `incident` must be preallocated - use mesh_count_incident().
// Returns a subset of a connectivity - entities may be repeated!
int32 mesh_get_incident(Mesh *mesh,
                        MeshConnectivity *incident, int32 dim,
                        Indices *entities, int32 dent);

// Get local ids of non-unique entities `incident` of dimension `dim` incident
// to `entities` of dimension `dent`, see `mesh_get_incident()`, with respect
// to `entities`. `local_ids` must be preallocated to same size as `incident`.
int32 mesh_get_local_ids(Mesh *mesh, Indices *local_ids,
                         Indices *entities, int32 dent,
                         MeshConnectivity *incident, int32 dim);

// Select entities of dimension `dim` that are completely given by entities of
// dimension `dent`.
// Example: mesh_select_complete(mesh, mask, 2, vertices, 0) will select all
// complete faces, whose vertices are listed in `vertices`.
int32 mesh_select_complete(Mesh *mesh, Mask *mask, int32 dim,
                           Indices *entities, int32 dent);

// Return the coordinates of centroids of mesh entities with dimension `dim`.
// `ccoors` must be preallocated.
int32 mesh_get_centroids(Mesh *mesh, float64 *ccoors, int32 dim);

// Return the volumes of mesh entities with dimension `dim` > 0.
int32 mesh_get_volumes(Mesh *mesh, float64 *volumes, int32 dim);

// Return the normals of facets for each mesh cell. The normals can be accessed
// using the cell-facet connectivity.
// If `which` is -1, two normals of each quadrilateral face are averaged. If it
// is 0 or 1, the corresponding normal is used.
int32 mesh_get_facet_normals(Mesh *mesh, float64 *normals, int32 which);

int32 me_get_incident(MeshEntity *entity, Indices *out, int32 dim);
int32 me_get_incident2(MeshEntity *entity, Indices *out,
                       MeshConnectivity *conn);
int32 contains(Indices *i1, Indices *i2);

int32 get_local_connectivity(MeshConnectivity *loc,
                             Indices *cell_vertices,
                             MeshConnectivity *refloc);
int32 sort_local_connectivity(MeshConnectivity *loc, uint32 *oris,
                              uint32 num);
void uint32_sort234_copy(uint32 *out, uint32 *p, uint32 num);
int32 uint32_sort4(uint32 *p);
int32 uint32_sort3(uint32 *p);
int32 uint32_sort2(uint32 *p);

// Mesh entity is given by (dimension, index) or (dim, ii).
// Index ii is in [0, num[dim] - 1].
int32 mesh_get_entity(int32 *n_items, int32 dim, int32 ii);

// ?
int32 mesh_iterate_entities(Mesh *mesh, int32 *entity, int32 dim, Region *reg);

#endif /* !MESH_H */

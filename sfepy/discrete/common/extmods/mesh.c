#include "common.h"
#include "geomtrans.h"
#include "mesh.h"

static void debprintf(const char *what, ...)
{
#ifdef DEBUG_MESH
  va_list ap;

  va_start(ap, what);
  vprintf(what, ap);
  va_end(ap);
#endif
}

int32 mesh_init(Mesh *mesh)
{
  uint32 ii;
  MeshTopology *topology = 0;
  LocalEntities *entities = 0;

  topology = mesh->topology;

  topology->max_dim = 0;
  topology->cell_types = 0;
  topology->edge_oris = 0;
  topology->face_oris = 0;
  memset(topology->num, 0, 16 * sizeof(int32));
  memset(topology->_conn, 0, 16 * sizeof(MeshConnectivity));
  for (ii = 0; ii < 16; ii++) {
    topology->conn[ii] = &topology->_conn[ii];
    topology->conn[ii]->num = 0;
    topology->conn[ii]->indices = 0;
    topology->conn[ii]->offsets = 0;
  }

  mesh->geometry->num = 0;
  mesh->geometry->dim = 0;
  mesh->geometry->coors = 0;

  entities = mesh->entities;
  entities->num = MAX_EL_TYPES;
  memset(entities->_edges, 0, MAX_EL_TYPES * sizeof(MeshConnectivity));
  for (ii = 0; ii < MAX_EL_TYPES; ii++) {
    entities->edges[ii] = &entities->_edges[ii];
    entities->edges[ii]->num = 0;
    entities->edges[ii]->indices = 0;
    entities->edges[ii]->offsets = 0;
  }
  memset(entities->_faces, 0, MAX_EL_TYPES * sizeof(MeshConnectivity));
  for (ii = 0; ii < MAX_EL_TYPES; ii++) {
    entities->faces[ii] = &entities->_faces[ii];
    entities->faces[ii]->num = 0;
    entities->faces[ii]->indices = 0;
    entities->faces[ii]->offsets = 0;
  }

  return(RET_OK);
}

#undef __FUNC__
#define __FUNC__ "mesh_free"
int32 mesh_free(Mesh *mesh)
{
  int32 ii;
  MeshTopology *topology = mesh->topology;
  LocalEntities *entities = mesh->entities;

  for (ii = 0; ii < 16; ii++) {
    conn_free(topology->conn[ii]);
  }

  for (ii = 0; ii < MAX_EL_TYPES; ii++) {
    conn_free(entities->edges[ii]);
    conn_free(entities->faces[ii]);
  }

  free_mem(topology->edge_oris);
  free_mem(topology->face_oris);

  return(RET_OK);
}

int32 mesh_print(Mesh *mesh, FILE *file, int32 header_only)
{
  uint32 ii, id;
  uint32 *num;
  uint32 D = mesh->topology->max_dim;
  MeshGeometry *geometry = mesh->geometry;
  MeshTopology *topology = mesh->topology;
  MeshConnectivity *conn = 0;

  fprintf(file, "Mesh %p (vertices: %d dimension: %d)\n",
          mesh, geometry->num, geometry->dim);
  fprintf(file, "topology: max_dim: %d\n", topology->max_dim);
  num = topology->num;
  fprintf(file, "n_cell: %d, n_face: %d, n_edge: %d, n_vertex: %d\n",
          num[3], num[2], num[1], num[0]);

  if (header_only == 0) { // Full print.
    fprintf(file, "vertex coordinates:\n");
    for (ii = 0; ii < geometry->num; ii++) {
      for (id = 0; id < geometry->dim; id++) {
        fprintf(file, " %.8e", geometry->coors[geometry->dim * ii + id]);
      }
      fprintf(file, "\n");
    }

    fprintf(file, "topology connectivities:\n");
    for (ii = 0; ii <= D; ii++) {
      for (id = 0; id <= D; id++) {
        fprintf(file, "incidence %d -> %d:\n", ii, id);
        conn = topology->conn[IJ(D, ii, id)];
        conn_print(conn, file);
      }
    }
  }

  return(RET_OK);
}

inline int32 mei_init(MeshEntityIterator *iter, Mesh *mesh, uint32 dim)
{
  iter->entity->mesh = mesh;
  iter->entity->dim = dim;
  iter->entity->ii = 0;
  iter->it = 0;
  iter->ptr = 0;
  iter->it_end = mesh->topology->num[dim];

  return(RET_OK);
}

inline int32 mei_init_conn(MeshEntityIterator *iter, MeshEntity *entity,
                           uint32 dim)
{
  Mesh *mesh = entity->mesh;
  uint32 D = mesh->topology->max_dim;
  MeshConnectivity *conn = mesh->topology->conn[IJ(D, entity->dim, dim)];

  iter->entity->mesh = mesh;
  iter->entity->dim = dim;
  iter->it = 0;

  if ((conn->num > 0) && (conn->indices > 0)) {
    iter->ptr = conn->indices + conn->offsets[entity->ii];
    iter->it_end = conn->offsets[entity->ii+1] - conn->offsets[entity->ii];
    iter->entity->ii = iter->ptr[iter->it];
  } else {
    iter->ptr = 0;
    iter->it_end = 0;
    iter->entity->ii = 0;
  }

  return(RET_OK);
}

inline int32 mei_init_sub(MeshEntityIterator *iter, Mesh *mesh,
                          Indices *entities, uint32 dim)
{
  iter->entity->mesh = mesh;
  iter->entity->dim = dim;
  iter->it = 0;

  iter->ptr = entities->indices;
  iter->it_end = entities->num;
  iter->entity->ii = iter->ptr[iter->it];

  return(RET_OK);
}

int32 mei_print(MeshEntityIterator *iter, FILE *file)
{
  fprintf(file, "it: %d, entity: dim: %d, ii: %d\n",
          iter->it, iter->entity->dim, iter->entity->ii);
  return(RET_OK);
}

inline int32 mei_go(MeshEntityIterator *iter)
{
  return(iter->it < iter->it_end);
}

inline int32 mei_next(MeshEntityIterator *iter)
{
  iter->it += 1;
  if (iter->it < iter->it_end) {
    iter->entity->ii = (iter->ptr ? iter->ptr[iter->it] : iter->it);
  }

  return(RET_OK);
}

int32 ind_print(Indices *ind, FILE *file)
{
  uint32 ii;
  if (!ind) return(RET_OK);

  fprintf(file, "indices: num: %d\n", ind->num);
  for (ii = 0; ii < ind->num; ii++) {
    fprintf(file, "%d: %d\n", ii, ind->indices[ii]);
  }

  return(RET_OK);
}

#undef __FUNC__
#define __FUNC__ "conn_alloc"
int32 conn_alloc(MeshConnectivity *conn, uint32 num, uint32 n_incident)
{
  int32 ret = RET_OK;

  if ((conn->num > 0) && (conn->n_incident > 0)) {
    conn_free(conn);
  }

  if (num > 0) {
    conn->num = num;
    conn->offsets = alloc_mem(uint32, num + 1);
    ERR_CheckGo(ret);
  }

  if (n_incident > 0) {
    conn->n_incident = n_incident;
    conn->indices = alloc_mem(uint32, n_incident);
    ERR_CheckGo(ret);
  } else if (num == 0) { // Empty connectivity.
    conn->n_incident = n_incident;
    conn->indices = 0;
  }

 end_label:
  if (ERR_Chk) {
    conn_free(conn);
  }

  return(ret);
}

#undef __FUNC__
#define __FUNC__ "conn_resize"
int32 conn_resize(MeshConnectivity *conn, uint32 num, uint32 n_incident)
{
  int32 ret = RET_OK;

  conn->num = num;
  conn->offsets = realloc_mem(conn->offsets, uint32, num + 1);
  ERR_CheckGo(ret);

  conn->n_incident = n_incident;
  conn->indices = realloc_mem(conn->indices, uint32, n_incident);
  ERR_CheckGo(ret);

 end_label:
  if (ERR_Chk) {
    errput("conn_resize() failed!");
  }

  return(ret);
}

#undef __FUNC__
#define __FUNC__ "conn_free"
int32 conn_free(MeshConnectivity *conn)
{
   free_mem(conn->indices);
   free_mem(conn->offsets);
   conn->num = 0;
   conn->n_incident = 0;

   return(RET_OK);
}

int32 conn_print(MeshConnectivity *conn, FILE *file)
{
  uint32 ii, ic;
  if (!conn) return(RET_OK);

  fprintf(file, "conn: num: %d, n_incident: %d\n", conn->num, conn->n_incident);

  for (ii = 0; ii < conn->num; ii++) {
    fprintf(file, "%d:", ii);
    for (ic = conn->offsets[ii]; ic < conn->offsets[ii+1]; ic++) {
      fprintf(file, " %d", conn->indices[ic]);
    }
    fprintf(file, "\n");
  }

  return(RET_OK);
}

// Set offsets and indices of conn from other.
inline int32 conn_set_from(MeshConnectivity *conn, MeshConnectivity *other)
{
  memcpy(conn->offsets, other->offsets, (conn->num + 1) * sizeof(uint32));
  memcpy(conn->indices, other->indices, conn->n_incident * sizeof(uint32));

  return(RET_OK);
}

inline int32 conn_set_to_free(MeshConnectivity *conn, uint32 ii,
                              uint32 incident)
{
  uint32 ok;
  uint32 *off, *ptr;

  off = conn->offsets + ii;
  ptr = conn->indices + off[0];
  ok = 0;
  while (ptr < (conn->indices + off[1])) {
    if (ptr[0] == UINT32_None) { // Not found & free slot.
      ptr[0] = incident;
      ok = 1;
      break;
    }
    ptr++;
  }
  if (!ok) {
    errput("no free connectivity position (internal error)!\n");
    return(RET_Fail);
  } else {
    return(RET_OK);
  }
}

int32 mesh_set_coors(Mesh *mesh, float64 *coors, int32 num, int32 dim,
                     int32 tdim)
{
  MeshGeometry *geometry = mesh->geometry;

  geometry->coors = coors;
  geometry->num = num;
  geometry->dim = dim;

  mesh->topology->max_dim = tdim;
  mesh->topology->num[0] = num;

  return(RET_OK);
}

int32 mesh_setup_connectivity(Mesh *mesh, int32 d1, int32 d2)
{
  int32 ret = RET_OK;
  int32 d3 = 0;
  MeshTopology *topology = mesh->topology;
  uint32 D = topology->max_dim;

  debprintf("request connectivity %d -> %d\n", d1, d2);

  if (topology->num[d1] == 0) {
    mesh_build(mesh, d1);
    ERR_CheckGo(ret);
  }

  if (topology->num[d2] == 0) {
    mesh_build(mesh, d2);
    ERR_CheckGo(ret);
  }

  if (topology->conn[IJ(D, d1, d2)]->num) {
    return(ret);
  }

  if (d1 < d2) {
    mesh_setup_connectivity(mesh, d2, d1);
    mesh_transpose(mesh, d1, d2);
  } else {
    if ((d1 == 0) && (d2 == 0)) {
      d3 = D;
    } else {
      d3 = 0;
    }
    if ((d1 > 0) && (d2 == 0)) {
      errput("connectivity %d -> %d should already exist!\n", d1, d2);
      ERR_CheckGo(ret);
    }
    mesh_setup_connectivity(mesh, d1, d3);
    mesh_setup_connectivity(mesh, d3, d2);
    mesh_intersect(mesh, d1, d2, d3);
  }
  ERR_CheckGo(ret);

 end_label:
  return(ret);
}

int32 mesh_free_connectivity(Mesh *mesh, int32 d1, int32 d2)
{
  uint32 D = mesh->topology->max_dim;
  MeshConnectivity *conn = 0;

  debprintf("free connectivity %d -> %d\n", d1, d2);

  conn = mesh->topology->conn[IJ(D, d1, d2)];
  conn_free(conn);

  return(RET_OK);
}

int32 mesh_build(Mesh *mesh, int32 dim)
{
  int32 ret = RET_OK;
  uint32 n_incident, n_v_max, n_loc;
  uint32 ii, ic, id, found;
  uint32 D = mesh->topology->max_dim;
  uint32 facet[4]; // Max. space for single facet.
  uint32 *oris = 0;
  uint32 loc_oris[12];
  uint32 *nDd = 0;
  uint32 *ptr1 = 0, *ptr2 = 0;
  uint32 *cell_types = mesh->topology->cell_types;
  Indices cell_vertices[1];
  MeshEntityIterator it0[1], it1[1], it2[1];
  MeshConnectivity *cD0 = 0; // D -> 0 - known
  MeshConnectivity *cDd = 0; // D -> d - to compute
  MeshConnectivity *cd0 = 0; // d -> 0 - to compute
  MeshConnectivity **locs = 0;
  MeshConnectivity *loc = 0;
  MeshConnectivity sloc[1]; // Local connectivity with sorted global vertices.
  MeshConnectivity gloc[1]; // Local connectivity with global vertices.

  debprintf("build %d\n", dim);

  if (!mesh->topology->conn[IJ(D, D, D)]->num) {
    mesh_setup_connectivity(mesh, D, D);
  }

  cD0 = mesh->topology->conn[IJ(D, D, 0)];
  cDd = mesh->topology->conn[IJ(D, D, dim)];
  cd0 = mesh->topology->conn[IJ(D, dim, 0)];

  locs = (dim == 1) ? mesh->entities->edges : mesh->entities->faces;

  // Max. number of vertices in facets.
  n_v_max = (dim == 1) ? 2 : 4;

  // Count entities of D -> d.
  conn_alloc(cDd, mesh->topology->num[D], 0);
  ERR_CheckGo(ret);
  nDd = cDd->offsets + 1;

  for (mei_init(it0, mesh, D); mei_go(it0); mei_next(it0)) {
    loc = locs[cell_types[it0->it]];
    nDd[it0->it] = loc->num;
  }

  // cDd->offsets now contains counts - make a cumsum to get offsets.
  for (ii = 1; ii < cDd->num + 1; ii++) {
    cDd->offsets[ii] += cDd->offsets[ii-1];
  }

  n_incident = cDd->offsets[cDd->num];
  debprintf("build n_incident (%d -> %d): %d\n", D, dim, n_incident);

  // Cell-local orientations w.r.t. D -> d.
  oris = alloc_mem(uint32, n_incident);
  if (dim == 2) {
    free_mem(mesh->topology->face_oris);
    mesh->topology->face_oris = oris;
  } else {
    free_mem(mesh->topology->edge_oris);
    mesh->topology->edge_oris = oris;
  }

  // Allocate D -> d indices.
  conn_alloc(cDd, 0, n_incident);
  ERR_CheckGo(ret);
  for (ii = 0; ii < cDd->n_incident; ii++) {
    cDd->indices[ii] = UINT32_None; // "not set" value.
  }

  // Allocate maximal buffers for d -> 0 arrays.
  conn_alloc(cd0, n_incident, n_incident * n_v_max);
  debprintf("build max. n_incident_vertex: %d\n", n_incident * n_v_max);

  // Allocate maximal buffers for local connectivity with sorted global
  // vertices. Counts have to be set to zero to avoid spurious calls to
  // conn_free()!
  sloc->num = sloc->n_incident = 0;
  conn_alloc(sloc, 12, n_v_max * 12);
  gloc->num = gloc->n_incident = 0;
  conn_alloc(gloc, 12, n_v_max * 12);

  id = 0;
  for (mei_init(it0, mesh, D); mei_go(it0); mei_next(it0)) {
    // Get vertex sets for local entities of current cell.
    loc = locs[cell_types[it0->it]];
    me_get_incident2(it0->entity, cell_vertices, cD0);
    get_local_connectivity(gloc, cell_vertices, loc);
    conn_set_from(sloc, gloc);
    sort_local_connectivity(sloc, loc_oris, loc->num);

    // Iterate over entities in the vertex sets.
    for (ii = 0; ii < loc->num; ii++) {
      // ii points to a vertex set in sloc/gloc.
      n_loc = sloc->offsets[ii+1] - sloc->offsets[ii];

      // Try to find entity in cells visited previously.
      for (mei_init_conn(it1, it0->entity, D); mei_go(it1); mei_next(it1)) {
        if (it1->entity->ii >= it0->entity->ii) continue;

        // Iterate over facets of visited cells.
        for (mei_init_conn(it2, it1->entity, dim); mei_go(it2); mei_next(it2)) {
          ptr1 = cd0->indices + cd0->offsets[it2->entity->ii];
          uint32_sort234_copy(facet, ptr1, n_loc);
          ptr2 = sloc->indices + sloc->offsets[ii];
          found = 1;
          for (ic = 0; ic < n_loc; ic++) {
            if (facet[ic] != ptr2[ic]) {
              found = 0;
              break;
            }
          }
          if (found) {
            // Assign existing entity to D -> d.
            conn_set_to_free(cDd, it0->entity->ii, it2->entity->ii);
            goto found_label;
          }
        }
      }
      // Entity not found - create new.

      // Add it as 'id' to D -> d.
      conn_set_to_free(cDd, it0->entity->ii, id);

      // Add vertices in gloc to d -> 0.
      cd0->offsets[id+1] = cd0->offsets[id] + n_loc;

      ptr1 = cd0->indices + cd0->offsets[id];
      ptr2 = gloc->indices + gloc->offsets[ii];
      for (ic = 0; ic < n_loc; ic++) {
        ptr1[ic] = ptr2[ic];
      }

      // Increment entity counter.
      id++;

    found_label:
      // Store entity orientation key to position of the last used item in cDd.
      ptr1 = cDd->offsets + it0->entity->ii;
      ic = ptr1[1] - 1;
      while (ic >= ptr1[0]) {
        if (cDd->indices[ic] != UINT32_None) { // Not found & free slot.
          break;
        }
        ic--;
      }
      // printf("%d << %d, %d\n", ic, ii, loc_oris[ii]);
      oris[ic] = loc_oris[ii];
    }
  }
  debprintf("build n_unique: %d, n_incident (%d -> 0): %d\n",
            id, dim, cd0->offsets[id]);

  // Update entity count in topology.
  mesh->topology->num[dim] = id;

  // Strip d -> 0.
  conn_resize(cd0, id, cd0->offsets[id]);

 end_label:
  conn_free(sloc);
  conn_free(gloc);

  return(ret);
}

#undef __FUNC__
#define __FUNC__ "mesh_transpose"
int32 mesh_transpose(Mesh *mesh, int32 d1, int32 d2)
{
  int32 ret = RET_OK;
  uint32 n_incident;
  uint32 ii;
  uint32 *nd2 = 0;
  uint32 D = mesh->topology->max_dim;
  MeshEntityIterator it2[1], it1[1];
  MeshConnectivity *c12 = 0; // d1 -> d2 - to compute

  debprintf("transpose %d -> %d\n", d1, d2);

  if (d1 >= d2) {
    errput("d1 must be smaller than d2 in mesh_transpose()!\n");
    ERR_CheckGo(ret);
  }

  c12 = mesh->topology->conn[IJ(D, d1, d2)];

  // Count entities of d2 -> d1.
  conn_alloc(c12, mesh->topology->num[d1], 0);
  ERR_CheckGo(ret);
  nd2 = c12->offsets + 1;

  for (mei_init(it2, mesh, d2); mei_go(it2); mei_next(it2)) {
    for (mei_init_conn(it1, it2->entity, d1); mei_go(it1); mei_next(it1)) {
      nd2[it1->entity->ii]++;
    }
  }

  // c12->offsets now contains counts - make a cumsum to get offsets.
  for (ii = 1; ii < c12->num + 1; ii++) {
    c12->offsets[ii] += c12->offsets[ii-1];
  }

  n_incident = c12->offsets[c12->num];
  debprintf("transpose n_incident (%d -> %d): %d\n", d1, d2, n_incident);

  // Fill in the indices.
  conn_alloc(c12, 0, n_incident);
  ERR_CheckGo(ret);
  for (ii = 0; ii < c12->n_incident; ii++) {
    c12->indices[ii] = UINT32_None; // "not set" value.
  }

  for (mei_init(it2, mesh, d2); mei_go(it2); mei_next(it2)) {
    for (mei_init_conn(it1, it2->entity, d1); mei_go(it1); mei_next(it1)) {
      conn_set_to_free(c12, it1->entity->ii, it2->entity->ii);
      ERR_CheckGo(ret);
    }
  }

 end_label:
  return(ret);
}

#undef __FUNC__
#define __FUNC__ "mesh_intersect"
int32 mesh_intersect(Mesh *mesh, int32 d1, int32 d2, int32 d3)
{
  int32 ret = RET_OK;
  uint32 D = mesh->topology->max_dim;
  uint32 n_incident, ii;
  uint32 *nd2 = 0;
  char *mask = 0;
  MeshEntityIterator it1[1], it2[1], it3[1];
  Indices ei1[1], ei2[1];
  MeshConnectivity *c12 = 0; // d1 -> d2 - to compute
  MeshConnectivity *c10 = 0; // d1 -> 0 - known
  MeshConnectivity *c20 = 0; // d2 -> 0 - known

  debprintf("intersect %d -> %d (%d)\n", d1, d2, d3);

  if (d1 < d2) {
    errput("d1 must be greater or equal to d2 in mesh_intersect()!\n");
    ERR_CheckGo(ret);
  }

  c12 = mesh->topology->conn[IJ(D, d1, d2)];
  if (d1 > d2) {
    c10 = mesh->topology->conn[IJ(D, d1, 0)];
    c20 = mesh->topology->conn[IJ(D, d2, 0)];
  }

  mask = alloc_mem(char, mesh->topology->num[d2]);

  // Count entities of d2 -> d1.
  conn_alloc(c12, mesh->topology->num[d1], 0);
  ERR_CheckGo(ret);
  nd2 = c12->offsets + 1;

  for (mei_init(it1, mesh, d1); mei_go(it1); mei_next(it1)) {
    // Clear mask for it1 incident entities of dimension d2.
    for (mei_init_conn(it3, it1->entity, d3); mei_go(it3); mei_next(it3)) {
      for (mei_init_conn(it2, it3->entity, d2); mei_go(it2); mei_next(it2)) {
        mask[it2->entity->ii] = 0;
      }
    }

    for (mei_init_conn(it3, it1->entity, d3); mei_go(it3); mei_next(it3)) {
      for (mei_init_conn(it2, it3->entity, d2); mei_go(it2); mei_next(it2)) {
        if (mask[it2->entity->ii]) continue;
        mask[it2->entity->ii] = 1;

        if (d1 == d2) {
          if (it1->entity->ii != it2->entity->ii) {
            nd2[it1->entity->ii]++;
          }
        } else {
          // Get incident vertices.
          me_get_incident2(it1->entity, ei1, c10);
          me_get_incident2(it2->entity, ei2, c20);
          if (contains(ei1, ei2)) {
            nd2[it1->entity->ii]++;
          }
        }
      }
    }
  }

  // c12->offsets now contains counts - make a cumsum to get offsets.
  for (ii = 1; ii < c12->num + 1; ii++) {
    c12->offsets[ii] += c12->offsets[ii-1];
  }

  n_incident = c12->offsets[c12->num];
  debprintf("intersect n_incident (%d -> %d): %d\n", d1, d2, n_incident);

  // Fill in the indices.
  conn_alloc(c12, 0, n_incident);
  ERR_CheckGo(ret);
  for (ii = 0; ii < c12->n_incident; ii++) {
    c12->indices[ii] = UINT32_None; // "not set" value.
  }

  for (mei_init(it1, mesh, d1); mei_go(it1); mei_next(it1)) {
    // Clear mask for it1 incident entities of dimension d2.
    for (mei_init_conn(it3, it1->entity, d3); mei_go(it3); mei_next(it3)) {
      for (mei_init_conn(it2, it3->entity, d2); mei_go(it2); mei_next(it2)) {
        mask[it2->entity->ii] = 0;
      }
    }

    for (mei_init_conn(it3, it1->entity, d3); mei_go(it3); mei_next(it3)) {
      for (mei_init_conn(it2, it3->entity, d2); mei_go(it2); mei_next(it2)) {
        if (mask[it2->entity->ii]) continue;
        mask[it2->entity->ii] = 1;

        if (d1 == d2) {
          if (it1->entity->ii != it2->entity->ii) {
            conn_set_to_free(c12, it1->entity->ii, it2->entity->ii);
            ERR_CheckGo(ret);
          }
        } else {
          // Get incident vertices.
          me_get_incident2(it1->entity, ei1, c10);
          me_get_incident2(it2->entity, ei2, c20);
          if (contains(ei1, ei2)) {
            conn_set_to_free(c12, it1->entity->ii, it2->entity->ii);
            ERR_CheckGo(ret);
          }
        }
      }
    }
  }

 end_label:
  free_mem(mask);

  return(ret);
}

uint32 mesh_count_incident(Mesh *mesh, int32 dim,
                           Indices *entities, int32 dent)
{
  uint32 ii, num = 0;
  uint32 *ptr;
  uint32 D = mesh->topology->max_dim;
  MeshConnectivity *conn = mesh->topology->conn[IJ(D, dent, dim)];

  if (!conn->num) {
    errput("connectivity %d -> %d is not avaliable!\n", dent, dim);
    ERR_CheckGo(num);
  }

  for (ii = 0; ii < entities->num; ii++) {
    ptr = conn->offsets + entities->indices[ii];
    num += ptr[1] - ptr[0];
  }

 end_label:
  return(num);
}

// `incident` must be preallocated - use mesh_count_incident().
int32 mesh_get_incident(Mesh *mesh,
                        MeshConnectivity *incident, int32 dim,
                        Indices *entities, int32 dent)
{
  int32 ret = RET_OK;
  uint32 ii;
  uint32 D = mesh->topology->max_dim;
  MeshEntityIterator it0[1], it1[1];
  MeshConnectivity *conn = mesh->topology->conn[IJ(D, dent, dim)];

  if (!conn->num) {
    errput("connectivity %d -> %d is not avaliable!\n", dent, dim);
    ERR_CheckGo(ret);
  }

  ii = 0;
  incident->offsets[0] = 0;
  for (mei_init_sub(it0, mesh, entities, dent); mei_go(it0); mei_next(it0)) {
    for (mei_init_conn(it1, it0->entity, dim); mei_go(it1); mei_next(it1)) {
      incident->indices[ii++] = it1->entity->ii;
    }
    incident->offsets[it0->it + 1] = incident->offsets[it0->it] + it1->it_end;
  }

 end_label:
  return(ret);
}

// `local_ids` must be preallocated to same size as `entities`.
int32 mesh_get_local_ids(Mesh *mesh, Indices *local_ids,
                         Indices *entities, int32 dent,
                         MeshConnectivity *incident, int32 dim)
{
  int32 ret = RET_OK;
  uint32 ii, iind, ic, found;
  uint32 D = mesh->topology->max_dim;
  MeshEntity entity[1];
  MeshEntityIterator it1[1];
  MeshConnectivity *conn = mesh->topology->conn[IJ(D, dim, dent)];

  if (!conn->num) {
    errput("connectivity %d -> %d is not avaliable!\n", dim, dent);
    ERR_CheckGo(ret);
  }

  entity->mesh = mesh;
  entity->dim = dim;

  ii = 0;
  for (iind = 0; iind < incident->num; iind++) {
    for (ic = incident->offsets[iind]; ic < incident->offsets[iind+1]; ic++) {
      entity->ii = incident->indices[ic];
      // printf("%d: ? %d in %d\n", iind, entities->indices[iind], entity->ii);
      found = 0;
      for (mei_init_conn(it1, entity, dent); mei_go(it1); mei_next(it1)) {
        if (entities->indices[iind] == it1->entity->ii) {
          local_ids->indices[ii++] = it1->it;
          // printf("%d -> %d\n", ii, it1->it);
          found = 1;
          break; // Degenerate cases - 1. occurrence is returned.
        }
      }
      if (!found) {
        errput("entity (%d, %d) not found in entity (%d, %d)!\n",
               entities->indices[iind], dent, entity->ii, dim);
        ERR_CheckGo(ret);
      }
    }
  }

 end_label:
  return(ret);
}

#undef __FUNC__
#define __FUNC__ "mesh_select_complete"
// Allocates mask->mask.
int32 mesh_select_complete(Mesh *mesh, Mask *mask, int32 dim,
                           Indices *entities, int32 dent)
{
  int32 ret = RET_OK;
  uint32 D = mesh->topology->max_dim;
  uint32 ii, inum;
  char *ent_mask = 0;
  MeshEntityIterator it0[1], it1[1];
  MeshConnectivity *conn = mesh->topology->conn[IJ(D, dim, dent)];

  if (!conn->num) {
    errput("connectivity %d -> %d is not avaliable!\n", dim, dent);
    ERR_CheckGo(ret);
  }

  mask->mask = alloc_mem(char, conn->num);
  mask->num = conn->num;
  mask->n_true = 0;

  ent_mask = alloc_mem(char, mesh->topology->num[dent]);

  for (ii = 0; ii < entities->num; ii++) {
    ent_mask[entities->indices[ii]] = 1;
  }

  for (mei_init(it0, mesh, dim); mei_go(it0); mei_next(it0)) {
    inum = 0;
    for (mei_init_conn(it1, it0->entity, dent); mei_go(it1); mei_next(it1)) {
      if (ent_mask[it1->entity->ii]) inum++;
    }
    // Check if all entities with dimension dent incident to entity it0 are set
    // in ent_mask.
    if (inum == it1->it_end) {
      mask->mask[it0->it] = 1;
      mask->n_true++;
    }
  }

 end_label:
  free_mem(ent_mask);

  return(ret);
}

// `ccoors` must be preallocated.
int32 mesh_get_centroids(Mesh *mesh, float64 *ccoors, int32 dim)
{
  uint32 nc = mesh->geometry->dim;
  uint32 id;
  float64 *ptr = ccoors;
  float64 *coors = mesh->geometry->coors;
  MeshEntityIterator it0[1], it1[1];

  for (mei_init(it0, mesh, dim); mei_go(it0); mei_next(it0)) {
    for (id = 0; id < nc; id++) {
      ptr[id] = 0.0;
    }
    for (mei_init_conn(it1, it0->entity, 0); mei_go(it1); mei_next(it1)) {
      for (id = 0; id < nc; id++) {
        ptr[id] += coors[nc * it1->entity->ii + id];
      }
    }
    for (id = 0; id < nc; id++) {
      ptr[id] /= (float64) it1->it_end;
    }
    ptr += nc;
  }

  return(RET_OK);
}

static inline float64 _det3x3(float64 j[9])
{
  return (j[0]*j[4]*j[8] + j[3]*j[7]*j[2] + j[1]*j[5]*j[6]
          - j[2]*j[4]*j[6] - j[5]*j[7]*j[0] - j[1]*j[3]*j[8]);
}

static inline float64 _tri_area(float64 *coors, uint32 *indices, uint32 nc)
{
#define VS(ic, id) (coors[nc*indices[ic] + id])

  uint32 id;
  float64 vv0;
  float64 v0[3], v1[3], ndir[3];

  if (nc == 2) { // 2D.
    v0[2] = 0.0;
    v1[2] = 0.0;
  }

  for (id = 0; id < nc; id++) {
    vv0 = VS(0, id);
    v0[id] = VS(1, id) - vv0;
    v1[id] = VS(2, id) - vv0;
  }
  gtr_cross_product(ndir, v0, v1);
  if (nc == 2) {
    return 0.5 * fabs(ndir[2]);
  } else {
    return 0.5 * sqrt(ndir[0] * ndir[0] + ndir[1] * ndir[1]
                      + ndir[2] * ndir[2]);
  }

#undef VS
}

static inline float64 _aux_hex(float64 *coors, uint32 *indices, uint32 nc,
                               int32 h, int32 i, int32 j, int32 k)
{
#define VS(ic, id) (coors[nc*indices[ic] + id])

  float64 mtx[9];

  mtx[0] = VS(i, 0) + VS(h, 0);
  mtx[1] = VS(i, 1) + VS(h, 1);
  mtx[2] = VS(i, 2) + VS(h, 2);
  mtx[3] = VS(j, 0) - VS(i, 0);
  mtx[4] = VS(j, 1) - VS(i, 1);
  mtx[5] = VS(j, 2) - VS(i, 2);
  mtx[6] = VS(k, 0) - VS(h, 0);
  mtx[7] = VS(k, 1) - VS(h, 1);
  mtx[8] = VS(k, 2) - VS(h, 2);

  return _det3x3(mtx);

#undef VS
}

// `volumes` must be preallocated.
int32 mesh_get_volumes(Mesh *mesh, float64 *volumes, int32 dim)
{
#define VS(ic, id) (coors[nc*entity_vertices->indices[ic] + id])

  int32 ret = RET_OK;
  uint32 D = mesh->topology->max_dim;
  uint32 nc = mesh->geometry->dim;
  uint32 id;
  uint32 indx2[6];
  float64 vol, aux, vv0, vv1, vv2, vv3;
  float64 *ptr = volumes;
  float64 *coors = mesh->geometry->coors;
  float64 v0[3], v1[3], v2[3], ndir[3];
  Indices entity_vertices[1];
  MeshEntityIterator it0[1];
  MeshConnectivity *cd0 = 0; // d -> 0

  if (!dim) {
    errput("vertices have no volume!\n");
    ERR_CheckGo(ret);
  }

  cd0 = mesh->topology->conn[IJ(D, dim, 0)];

  for (mei_init(it0, mesh, dim); mei_go(it0); mei_next(it0)) {
    me_get_incident2(it0->entity, entity_vertices, cd0);
    /* mei_print(it0, stdout); */
    /* ind_print(entity_vertices, stdout); */
    vol = 0;

    if (dim == 1) { // Edges.
      for (id = 0; id < nc; id++) {
        aux = VS(1, id) - VS(0, id);
        vol += aux * aux;
      }
      ptr[0] = sqrt(vol);

    } else if (entity_vertices->num == 3) { // Triangles.
      ptr[0] = _tri_area(coors, entity_vertices->indices, nc);

    } else if (nc == 2) { // Quadrilateral cells.
      ptr[0] = _tri_area(coors, entity_vertices->indices, nc);
      indx2[0] = entity_vertices->indices[2];
      indx2[1] = entity_vertices->indices[3];
      indx2[2] = entity_vertices->indices[0];
      ptr[0] += _tri_area(coors, indx2, nc);

    } else if (nc == 3) { // 3D.
      if (entity_vertices->num == 4) {
        if (dim == 2) { // Quadrilateral faces (approximate).
          indx2[0] = entity_vertices->indices[0];
          indx2[1] = entity_vertices->indices[1];
          indx2[2] = entity_vertices->indices[2];
          indx2[3] = entity_vertices->indices[3];
          indx2[4] = entity_vertices->indices[0];
          indx2[5] = entity_vertices->indices[1];

          aux = _tri_area(coors, indx2, nc);
          aux += _tri_area(coors, indx2 + 1, nc);
          aux += _tri_area(coors, indx2 + 2, nc);
          aux += _tri_area(coors, indx2 + 3, nc);
          ptr[0] = 0.5 * aux;

        } else { // Tetrahedral cells.
          for (id = 0; id < nc; id++) {
            vv0 = VS(0, id);
            vv1 = VS(1, id);
            vv2 = VS(2, id);
            vv3 = VS(3, id);

            v0[id] = vv1 - vv0;
            v1[id] = vv2 - vv0;
            v2[id] = vv3 - vv2;
          }
          gtr_cross_product(ndir, v0, v1);
          gtr_dot_v3(ptr, v2, ndir, 3);
          ptr[0] /= 6.0;
        }
      } else if (entity_vertices->num == 6) { // Wedge cells
          for (id = 0; id < nc; id++) {
            vv0 = VS(0, id);
            vv1 = VS(1, id);
            vv2 = VS(2, id);
            vv3 = VS(3, id);

            v0[id] = vv1 - vv0;
            v1[id] = vv2 - vv0;
            v2[id] = vv3 - vv2;
          }
          gtr_cross_product(ndir, v0, v1);
          gtr_dot_v3(ptr, v2, ndir, 3);
          ptr[0] /= 2.0;
      } else { // Hexahedral cells with trilinear interpolation.
        // See https://math.stackexchange.com/questions/1628540/what-is-the-enclosed-volume-of-an-irregular-cube-given-the-x-y-z-coordinates-of
        // Uses 0 1 3 2 4 5 7 6 ordering w.r.t. sfepy.
        aux = _aux_hex(coors, entity_vertices->indices, nc, 0, 1, 3, 2);
        aux -= _aux_hex(coors, entity_vertices->indices, nc, 4, 5, 7, 6);
        aux -= _aux_hex(coors, entity_vertices->indices, nc, 0, 1, 4, 5);
        aux += _aux_hex(coors, entity_vertices->indices, nc, 3, 2, 7, 6);
        aux += _aux_hex(coors, entity_vertices->indices, nc, 0, 3, 4, 7);
        aux -= _aux_hex(coors, entity_vertices->indices, nc, 1, 2, 5, 6);
        ptr[0] = aux / 12.0;

      }
    }

    ptr += 1;
  }

 end_label:
  return(ret);

#undef VS
}

// `normals` must be preallocated.
int32 mesh_get_facet_normals(Mesh *mesh, float64 *normals, int32 which)
{
#define VS(ic, id) (coors[nc*cell_vertices->indices[ik[ic]] + id])

  uint32 D = mesh->topology->max_dim;
  uint32 nc = mesh->geometry->dim;
  int32 dim = D - 1;
  uint32 ii, id, n_loc;
  uint32 *ik;
  uint32 *cell_types = mesh->topology->cell_types;
  float64 *coors = mesh->geometry->coors;
  float64 vv0, vv1, vv2, vv3, v0[3], v1[3], v2[3], v3[3], ndir[3], ndir1[3];
  Indices cell_vertices[1];
  MeshEntityIterator it0[1];
  MeshConnectivity *cD0 = 0; // D -> 0
  MeshConnectivity *cDd = 0; // D -> d
  MeshConnectivity *loc = 0;
  MeshConnectivity **locs = 0;

  cD0 = mesh->topology->conn[IJ(D, D, 0)];
  cDd = mesh->topology->conn[IJ(D, D, dim)];

  // Local entities - reference cell edges or faces.
  locs = (dim == 1) ? mesh->entities->edges : mesh->entities->faces;

  for (mei_init(it0, mesh, D); mei_go(it0); mei_next(it0)) {
    me_get_incident2(it0->entity, cell_vertices, cD0);
    loc = locs[cell_types[it0->it]];

    for (ii = 0; ii < loc->num; ii++) {
      ik = loc->indices + loc->offsets[ii]; // Points to local facet vertices.

      n_loc = loc->offsets[ii+1] - loc->offsets[ii];
      if (n_loc == 2) { // Edge normals.
        for (id = 0; id < nc; id++) {
          v0[id] = VS(1, id) - VS(0, id);
        }
        ndir[0] = v0[1];
        ndir[1] = -v0[0];

      } else if (n_loc == 3) { // Triangular face normals.
        for (id = 0; id < nc; id++) {
          vv0 = VS(0, id);
          v0[id] = VS(1, id) - vv0;
          v1[id] = VS(2, id) - vv0;
        }
        gtr_cross_product(ndir, v0, v1);

      } else if (n_loc == 4) { // Quadrilateral face normals.
        for (id = 0; id < nc; id++) {
          vv0 = VS(0, id);
          vv1 = VS(1, id);
          vv2 = VS(2, id);
          vv3 = VS(3, id);

          v0[id] = vv1 - vv0;
          v1[id] = vv3 - vv0;
          v2[id] = vv3 - vv2;
          v3[id] = vv1 - vv2;
        }

        if (which == 0) {
          gtr_cross_product(ndir, v0, v1);

        } else if (which == 1) {
          gtr_cross_product(ndir, v2, v3);

        } else {
          gtr_cross_product(ndir, v0, v1);
          gtr_cross_product(ndir1, v2, v3);
          for (id = 0; id < nc; id++) {
            ndir[id] += ndir1[id];
          }
        }
      }

      gtr_normalize_v3(ndir, ndir, nc, 0);
      for (id = 0; id < nc; id++) {
        normals[nc * (cDd->offsets[it0->it] + ii) + id] = ndir[id];
      }
    }
  }

  return(RET_OK);

#undef VS
}

inline int32 me_get_incident(MeshEntity *entity, Indices *out, int32 dim)
{
  int32 ret = RET_OK;
  Mesh *mesh = entity->mesh;
  uint32 D = mesh->topology->max_dim;
  MeshConnectivity *conn = mesh->topology->conn[IJ(D, entity->dim, dim)];

  if (!conn->num) {
    errput("required connectivity is not avaliable!\n");
    ERR_CheckGo(ret);
  }
  out->indices = conn->indices + conn->offsets[entity->ii];
  out->num = conn->offsets[entity->ii + 1] - conn->offsets[entity->ii];

 end_label:
  return(ret);
}

inline int32 me_get_incident2(MeshEntity *entity, Indices *out,
                              MeshConnectivity *conn)
{
  out->indices = conn->indices + conn->offsets[entity->ii];
  out->num = conn->offsets[entity->ii + 1] - conn->offsets[entity->ii];

  return(RET_OK);
}

inline int32 contains(Indices *i1, Indices *i2)
{
  // Check if all indices in i2 are contained in i1.
  // Brute force search is used.
  uint32 ii1, ii2, v2, ok;

  for (ii2 = 0; ii2 < i2->num; ii2++) {
    v2 = i2->indices[ii2];
    ok = 0;
    for (ii1 = 0; ii1 < i1->num; ii1++) {
      if (i1->indices[ii1] == v2) {
        ok = 1;
        break;
      }
    }
    if (!ok) return(0);
  }

  return(1);
}

inline int32 get_local_connectivity(MeshConnectivity *loc,
                                    Indices *cell_vertices,
                                    MeshConnectivity *refloc)
{
  uint32 ii, ic;

  for (ii = 0; ii < (refloc->num + 1); ii++) {
    loc->offsets[ii] = refloc->offsets[ii];
  }

  for (ii = 0; ii < refloc->num; ii++) {
    for (ic = refloc->offsets[ii]; ic < refloc->offsets[ii+1]; ic++) {
      loc->indices[ic] = cell_vertices->indices[refloc->indices[ic]];
    }
  }

  return(RET_OK);
}

inline int32 sort_local_connectivity(MeshConnectivity *loc, uint32 *oris,
                                     uint32 num)
{
  int32 key = -1;
  uint32 ii, n_item;

  if (num == 0) {
    num = loc->num;
  }

  for (ii = 0; ii < num; ii++) {
    n_item = loc->offsets[ii+1] - loc->offsets[ii];
    switch (n_item) {
    case 2:
      key = uint32_sort2(loc->indices + loc->offsets[ii]);
      break;
    case 3:
      key = uint32_sort3(loc->indices + loc->offsets[ii]);
      break;
    case 4:
      key = uint32_sort4(loc->indices + loc->offsets[ii]);
      break;
    }
    oris[ii] = key;
  }

  return(RET_OK);
}

#define SwapValues(a, b, work) do {\
  (work) = (a); (a) = (b); (b) = (work);\
} while (0)

#define SORT4(p, work) do {\
  if ((p)[0] > (p)[1]) SwapValues((p)[0], (p)[1], (work));\
  if ((p)[1] > (p)[2]) SwapValues((p)[1], (p)[2], (work));\
  if ((p)[2] > (p)[3]) SwapValues((p)[2], (p)[3], (work));\
  if ((p)[0] > (p)[1]) SwapValues((p)[0], (p)[1], (work));\
  if ((p)[1] > (p)[2]) SwapValues((p)[1], (p)[2], (work));\
  if ((p)[0] > (p)[1]) SwapValues((p)[0], (p)[1], (work));\
} while (0)

#define SORT3(p, work) do {\
  if ((p)[0] > (p)[1]) SwapValues((p)[0], (p)[1], (work));\
  if ((p)[1] > (p)[2]) SwapValues((p)[1], (p)[2], (work));\
  if ((p)[0] > (p)[1]) SwapValues((p)[0], (p)[1], (work));\
} while (0)

#define SORT2(p, work) do {\
  if ((p)[0] > (p)[1]) SwapValues((p)[0], (p)[1], (work));\
} while (0)

inline void uint32_sort234_copy(uint32 *out, uint32 *p, uint32 num)
{
  uint32 ii;
  uint32 work;

  for (ii = 0; ii < num; ii++) {
    out[ii] = p[ii];
  }
  switch (num) {
  case 2:
    SORT2(out, work);
    break;
  case 3:
    SORT3(out, work);
    break;
  case 4:
    SORT4(out, work);
    break;
  }
}

inline int32 uint32_sort4(uint32 *p)
{
  int32 key = 0;
  uint32 work;

  if (p[0] < p[1]) key += 1;
  if (p[0] < p[2]) key += 2;
  if (p[1] < p[2]) key += 4;
  if (p[0] < p[3]) key += 8;
  if (p[1] < p[3]) key += 16;
  if (p[2] < p[3]) key += 32;

  SORT4(p, work);

  return(key);
}

inline int32 uint32_sort3(uint32 *p)
{
  int32 key = 0;
  uint32 work;

  if (p[0] < p[1]) key += 1;
  if (p[0] < p[2]) key += 2;
  if (p[1] < p[2]) key += 4;

  SORT3(p, work);

  return(key);
}

inline int32 uint32_sort2(uint32 *p)
{
  int32 key = 0;
  uint32 work;

  if (p[0] < p[1]) key += 1;

  SORT2(p, work);

  return(key);
}

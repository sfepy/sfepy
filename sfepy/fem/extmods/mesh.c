#include "common.h"
#include "mesh.h"

int32 mesh_init(Mesh *mesh)
{
  int32 ii;
  MeshTopology *topology = 0;

  topology = mesh->topology;

  topology->max_dim = 0;
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

  return(RET_OK);
}

int32 mesh_print(Mesh *mesh, FILE *file, int32 header_only)
{
  int32 ii, id;
  uint32 *num;
  int32 D = mesh->topology->max_dim;
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

int32 conn_iter_init(ConnIter *iter, MeshConnectivity *conn)
{
  iter->conn = conn;
  iter->ii = 0;

  if (iter->conn->num > 0) {
    iter->num = iter->conn->offsets[1];
    iter->ptr = iter->conn->indices;
  } else {
    iter->num = 0;
    iter->ptr = 0;
  }

  return(RET_OK);
}

int32 conn_iter_print(ConnIter *iter, FILE *file)
{
  fprintf(file, "ii: %d, num: %d, ptr: %p\n", iter->ii, iter->num, iter->ptr);

  return(RET_OK);
}

int32 conn_iter_next(ConnIter *iter)
{
  iter->ii += 1;
  if (iter->ii < iter->conn->num) {
    iter->ptr += iter->num;

    uint32 *off = iter->conn->offsets + iter->ii;
    iter->num = off[1] - off[0];
  }

  return(iter->ii < iter->conn->num);
}

int32 conn_print(MeshConnectivity *conn, FILE *file)
{
  ConnIter iter[1];

  if (!conn) return(RET_OK);

  fprintf(file, "conn: num: %d, n_incident: %d\n", conn->num, conn->n_incident);

  conn_iter_init(iter, conn);

  while (conn_iter_next(iter)) {
    int32 ii;
    fprintf(file, "%d:", iter->ii);
    for (ii = 0; ii < iter->num; ii++) {
      fprintf(file, " %d", iter->ptr[ii]);
    }
    fprintf(file, "\n");
  }

  return(RET_OK);
}

int32 entity_iter_init(EntityIter *iter, int32 dim, MeshTopology *topology)
{
  iter->topology = topology;
  iter->ii = 0;
  iter->dim = dim;

  return(RET_OK);
}

int32 entity_iter_print(EntityIter *iter, FILE *file)
{
  fprintf(file, "dim: %d, ii: %d\n", iter->dim, iter->ii);

  return(RET_OK);
}

int32 entity_iter_next(EntityIter *iter)
{
  iter->ii += 1;
  return (iter->ii < iter->topology->num[iter->dim]);
}

int32 mesh_set_coors(Mesh *mesh, float64 *coors, int32 num, int32 dim)
{
  MeshGeometry *geometry = mesh->geometry;

  geometry->coors = coors;
  geometry->num = num;
  geometry->dim = dim;

  mesh->topology->max_dim = dim;
  mesh->topology->num[0] = num;

  return(RET_OK);
}

int32 mesh_setup_connectivity(Mesh *mesh, int32 d1, int32 d2)
{
  int32 ret = RET_OK;
  int32 d3 = 0;
  MeshConnectivity *conn = 0;
  MeshTopology *topology = mesh->topology;
  int32 D = topology->max_dim;

  if (topology->num[d1] == 0) {
    mesh_build(mesh, d1);
  }

  if (topology->num[d2] == 0) {
    mesh_build(mesh, d2);
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
    mesh_setup_connectivity(mesh, d1, d3);
    mesh_setup_connectivity(mesh, d3, d2);
    mesh_intersect(mesh, d1, d2, d3);
  }

 end_label:
  return(ret);
}

int32 mesh_build(Mesh *mesh, int32 dim)
{
  int32 ret = RET_OK;
  uint32 kk = 0;
  int32 D = mesh->topology->max_dim;
  MeshConnectivity *c1 = 0; // D -> d
  MeshConnectivity *c2 = 0; // d -> 0

  c1 = mesh->topology->conn[IJ(D, D, dim)];
  c2 = mesh->topology->conn[IJ(D, dim, 0)];


 end_label:
  return(ret);
}

int32 mesh_transpose(Mesh *mesh, int32 d1, int32 d2)
{
  int32 ret = RET_OK;
  uint32 kk = 0;
  int32 D = mesh->topology->max_dim;
  MeshConnectivity *c12 = 0; // d1 -> d2 - to compute
  MeshConnectivity *c21 = 0; // d2 -> d1 - known

  c12 = mesh->topology->conn[IJ(D, d1, d2)];
  c21 = mesh->topology->conn[IJ(D, d2, d1)];

 end_label:
  return(ret);
}

int32 mesh_intersect(Mesh *mesh, int32 d1, int32 d2, int32 d3)
{
  int32 ret = RET_OK;

 end_label:
  return(ret);
}

#include <stdlib.h>

#include "meshutils.h"
#include "sort.h"
#include "geomtrans.h"

int32 *ic;
int32 n_c;

int compare_row_cols( const void *obj_a, const void *obj_b )
{
  int32 *row_a, *row_b;
  int32 ii;
  int ret = 0;

  row_a = (int32 *) obj_a;
  row_b = (int32 *) obj_b;
  
/*   output( "%p %d\n", ic, n_c ); */
  for (ii = 0; ii < n_c; ii++) {
/*     output( "%d %d\n", row_a[ic[ii]], row_b[ic[ii]]  ); */
/*     output( "%d\n", ic[ii] ); */
    if (row_a[ic[ii]] < row_b[ic[ii]]) {
      ret = -1;
      break;
    } else if (row_a[ic[ii]] > row_b[ic[ii]]) {
      ret = 1;
      break;
    }
  }
  
  return( ret );
}

/*!
  @par Revision history:
  - 07.02.2006
*/
int32 sort_rows( int32 *array, int32 n_row, int32 n_col,
		int32 *i_sort_col, int32 n_sort_col )
{
/*   output( "%p %p\n", array, i_sort_col ); */
  if ((array == 0) || (i_sort_col == 0)) {
    output( "null pointer: %p %p\n", array, i_sort_col );
  }

  int32_sort_rows( array, n_row, n_col, i_sort_col, n_sort_col );
/*   ic = i_sort_col; */
/*   n_c = n_sort_col; */

/*   qsort( array, n_row, n_col * sizeof( int32 ), &compareRowCols ); */

/*   { */
/*     int32 ir, ic; */

/*     for (ir = 0; ir < n_row; ir++) { */
/*       for (ic = 0; ic < n_col; ic++) { */
/* 	output( "%d ", array[n_col*ir+ic] ); */
/*       } */
/*       output( "\n" ); */
/*     } */
/*   } */
/*   array[0] = -1; */

  return( RET_OK );
}

#define SwapValues( a, b, work ) do {\
  (work) = (a); (a) = (b); (b) = (work);\
} while (0)

#define MESH_Sort4( p, work ) do {\
  if ((p)[0] > (p)[1]) SwapValues( (p)[0], (p)[1], (work) );\
  if ((p)[1] > (p)[2]) SwapValues( (p)[1], (p)[2], (work) );\
  if ((p)[2] > (p)[3]) SwapValues( (p)[2], (p)[3], (work) );\
  if ((p)[0] > (p)[1]) SwapValues( (p)[0], (p)[1], (work) );\
  if ((p)[1] > (p)[2]) SwapValues( (p)[1], (p)[2], (work) );\
  if ((p)[0] > (p)[1]) SwapValues( (p)[0], (p)[1], (work) );\
} while (0)

#define MESH_Sort3( p, work ) do {\
  if ((p)[0] > (p)[1]) SwapValues( (p)[0], (p)[1], (work) );\
  if ((p)[1] > (p)[2]) SwapValues( (p)[1], (p)[2], (work) );\
  if ((p)[0] > (p)[1]) SwapValues( (p)[0], (p)[1], (work) );\
} while (0)

#define MESH_Sort2( p, work ) do {\
  if ((p)[0] > (p)[1]) SwapValues( (p)[0], (p)[1], (work) );\
} while (0)

/*!
  @par Revision history:
  - 14.01.2005, c
  - 17.01.2005
*/
int32 create_list( int32 *p_ii, int32 *objs, int32 objs_n_row, int32 objs_n_col,
		  int32 ig,
		  int32 *conn0, int32 conn_n_row, int32 conn_n_col,
		  int32 *items, int32 items_n_row, int32 items_n_col,
		  int32 is_sort )
{
  int32 *conn;
  int32 ii0, ii, ie, ij, io, iop, aux;

/*   output( "%d: %d %d %d %d %d %d\n", *p_ii, objs_n_row, objs_n_col, */
/* 	  conn_n_row, conn_n_col, items_n_row, items_n_col ); */
/*   output( "%p\n", conn0 ); */

  conn = conn0;
  ii0 = ii = *p_ii;
/*   output( "%d\n", ii ); */
  for (ie = 0; ie < conn_n_row; ie++) {
    for (ij = ii; ij < (ii + items_n_row); ij++) {
      io = ij - ii;
      objs[objs_n_col * ij + 0] = ig;
      objs[objs_n_col * ij + 1] = ie;
      objs[objs_n_col * ij + 2] = io;

      for (iop = 0; iop < items_n_col; iop++) {
/* 	output( "%d %d %d\n", ij, iop, conn[items[items_n_col*io+iop]] ); */
	aux = items[items_n_col*io+iop];
	if (aux >= 0)
	  objs[objs_n_col * ij + 3 + iop] = conn[aux];
      }
    }
    ii += items_n_row;
    conn += conn_n_col;
  }
  *p_ii = ii;
/*   output( "%d\n", ii ); */

  if (is_sort) {
    objs += ii0 * objs_n_col;
    if (items_n_col == 2) {
      for (ii = ii0; ii < (*p_ii); ii++) {
	MESH_Sort2( objs + 3, aux );
	objs += objs_n_col;
      }
    } else if (items_n_col == 3) {
      for (ii = ii0; ii < (*p_ii); ii++) {
	MESH_Sort3( objs + 3, aux );
	objs += objs_n_col;
      }
    } else if (items_n_col == 4) {
      for (ii = ii0; ii < (*p_ii); ii++) {
	if (objs[6] == -1) {
	  MESH_Sort3( objs + 3, aux );
	} else {
	  MESH_Sort4( objs + 3, aux );
	}
	objs += objs_n_col;
      }
    } else {
      return( RET_Fail );
    }
  }
/*   output( "%d\n", *p_ii ); */
  
  return( RET_OK );
}

/*!
  @par Revision history:
  - 18.01.2005, c
*/
int32 neighbour_list_ptr( int32 *pobj, int32 pobj_n_row,
			int32 *pg, int32 pg_n_row,
			int32 *pel, int32 pel_n_row,
			int32 *data_s, int32 data_s_n_row, int32 data_s_n_col,
			int32 *ic, int32 ic_n_row,
			int32 mode )
{
  int32 ii, i1, ip, icmp, cmp, n_cmp;
  int32 *aux0, *aux1;

  if (pobj_n_row == 0) return( RET_OK );

  n_cmp = (mode == 0) ? 2 : 4;

  pobj[0] = 0;
  ii = 0;
  while (ii < (data_s_n_row - 2)) {
    if (ic[ii]) {
      ip = ii;
      aux0 = data_s + data_s_n_col * ii + 3;
      while (1) {
	aux1 = data_s + data_s_n_col * (ip + 1) + 3;
	cmp = 1;
	for (icmp = 0; icmp < n_cmp; icmp++) {
	  if (aux0[icmp] != aux1[icmp]) {
	    cmp = 0;
	    break;
	  }
	}
	if (!cmp) break;
	ip += 1;
      }
/*       output( "neighbour_list_ptr %d %d\n", ii, ip ); */
      for (i1 = ii; i1 < (ip + 1); i1++) {
	aux1 = data_s + data_s_n_col * i1;
	pobj[pel[pg[aux1[0]] + aux1[1]] + aux1[2] + 1] += ip - ii;
      }
      ii = ip + 1;
    } else {
      ii += 1;
    }
  }
  return( RET_OK );
}

/*!
  @par Revision history:
  - 18.01.2005, c
  - 31.10.2005
*/
int32 neighbour_list( int32 *p_iu, int32 *objs, int32 objs_n_row,
		     int32 *uid, int32 uid_n_row,
		     int32 *cnt, int32 cnt_n_row,
		     int32 *pg, int32 pg_n_row,
		     int32 *pel, int32 pel_n_row,
		     int32 *pobj, int32 pobj_n_row,
		     int32 *data_s, int32 data_s_n_row, int32 data_s_n_col,
		     int32 *uid_in, int32 uid_in_n_row,
		     int32 *ic, int32 ic_n_row,
		     int32 *perm, int32 perm_n_row,
		     int32 mode )
{
  int32 ii, iu, i1, i2, ip, icmp, cmp, n_cmp, ptr;
  int32 *aux0, *aux1;

  if (objs_n_row == 0) return( RET_OK );

  n_cmp = (mode == 0) ? 2 : 4;

  ii = 0;
  iu = 0;
  while (ii < (data_s_n_row - 2)) {
    if (ic[ii]) {
      ip = ii;
      aux0 = data_s + data_s_n_col * ii + 3;
      while (1) {
	aux1 = data_s + data_s_n_col * (ip + 1) + 3;
	cmp = 1;
	for (icmp = 0; icmp < n_cmp; icmp++) {
	  if (aux0[icmp] != aux1[icmp]) {
	    cmp = 0;
	    break;
	  }
	}
	if (!cmp) break;
	ip += 1;
      }
/*       output( "neighbour_list %d %d\n", ii, ip ); */
      for (i1 = ii; i1 < (ip + 1); i1++) {
	// Own edge/face first.
	aux1 = data_s + data_s_n_col * i1;
	ptr = pel[pg[aux1[0]] + aux1[1]] + aux1[2];
	objs[pobj[ptr]+cnt[ptr]] = perm[i1];
	cnt[ptr] += 1;
	uid[ptr] = iu;
	if (iu != uid_in[i1]) {
	  errput( "%d == %d\n", iu, uid_in[i1] );
	  return( RET_Fail );
	}

	for (i2 = ii; i2 < (ip + 1); i2++) {
	  if (i1 == i2) continue;
	  objs[pobj[ptr]+cnt[ptr]] = perm[i2];
	  cnt[ptr] += 1;
	}
      }
      ii = ip + 1;
    } else {
      aux1 = data_s + data_s_n_col * ii;
      ptr = pel[pg[aux1[0]] + aux1[1]] + aux1[2];
      objs[pobj[ptr]+cnt[ptr]] = perm[ii];
      cnt[ptr] += 1;
      uid[ptr] = iu;
      if (iu != uid_in[ii]) {
	errput( "%d == %d\n", iu, uid_in[ii] );
	return( RET_Fail );
      }
      ii += 1;
    }
    iu += 1;
  }
  *p_iu = iu;

  return( RET_OK );
}

/*!
  No shape checks!

  @par Revision history:
  - 07.10.2005, c
*/
int32 orient_edges( int32 *ori, int32 ori_n_row, int32 ori_n_col,
		   int32 *conn, int32 conn_n_row, int32 conn_n_col,
		   int32 *edges, int32 edges_n_row, int32 edges_n_col )
{
  int32 ir, ic;

  for (ir = 0; ir < ori_n_row; ir++) {
    for (ic = 0; ic < ori_n_col; ic++) {
      ori[ic] = conn[edges[2*ic+1]] < conn[edges[2*ic]];
    }
    ori += ori_n_col;
    conn += conn_n_col;
  }

  return( RET_OK );
}

/*!
  No shape checks!

  @par Revision history:
  - 10.10.2005, c
*/
int32 assign_edge_nodes( int32 *p_iseq,  
		       int32 *econn, int32 econn_n_row, int32 econn_n_col,
		       int32 *cnt_en, int32 cnt_en_n_row, int32 cnt_en_n_col,
		       int32 *ori, int32 ori_n_row, int32 ori_n_col,
		       int32 *ntt, int32 ntt_n_row, int32 ntt_n_col,
		       int32 *uid, int32 uid_n_row,
		       int32 **edge_desc, int32 *edge_desc_n_row, int32 n_edge,
		       int32 cptr0 )
{
  int32 iel, ii, en, iep, it, ui, iseq, n_en;
  int32 *ei, *pori, *puid, *peconn, *pcnt;

  iseq = *p_iseq;
  pori = ori;
  puid = uid + cptr0;
  peconn = econn;
/*   output( "%d %d %d %d %d %d %d %d %d %d %d %d\n", */
/* 	  econn_n_row, econn_n_col, cnt_en_n_row, cnt_en_n_col, */
/* 	  ori_n_row, ori_n_col, ntt_n_row, ntt_n_col, */
/* 	  uid_n_row, n_edge, cptr0, iseq ); */
  for (iel = 0; iel < econn_n_row; iel++) {

    for (ii = 0; ii < n_edge; ii++) {
      ei = edge_desc[ii];
      n_en = edge_desc_n_row[ii];
      for (en = 0; en < n_en; en++) {
	iep = ei[2*en];
	it = ntt[ntt_n_col*pori[ii]+en];
	ui = puid[ii];
	pcnt = cnt_en + cnt_en_n_col * it + ui;
/* 	output( "%d %d %d %d %d %d %d\n", iel, ii, en, iep, it, ui, iseq ); */

	if ((*pcnt) >= 0) {
	  peconn[iep] = *pcnt;
	} else {
	  (*pcnt) = iseq;
	  peconn[iep] = iseq;
	  iseq++;
	}
      }
    }
    pori += ori_n_col;
    puid += n_edge;
    peconn += econn_n_col;
  }

  *p_iseq = iseq;

  return( RET_OK );
}

/*!
  No shape checks!

  @par Revision history:
  - 10.10.2005, c
  - 12.10.2005
*/
int32 interp_vertex_data( float64 *nod_out, int32 nod_out_n_row, int32 nod_out_n_col,
			int32 *econn, int32 econn_n_row, int32 econn_n_col,
			float64 *nod_in, int32 nod_in_n_row, int32 nod_in_n_col,
			int32 *conn, int32 conn_n_row, int32 conn_n_col,
			float64 *bf, int32 bf_n_row, int32 bf_n_col,
			int32 omit_cols )
{
  int32 iel, ir, ic, ii;
  int32 *pconn = conn, *peconn = econn;
  float64 val;

  for (iel = 0; iel < econn_n_row; iel++) {
    for (ir = 0; ir < bf_n_row; ir++){
      for (ic = 0; ic < (nod_out_n_col - omit_cols); ic++) {
	val = 0.0;
	for (ii = 0; ii < bf_n_col; ii++) {
	  val += bf[bf_n_col*ir+ii] * nod_in[nod_in_n_col*pconn[ii]+ic];
	}
	nod_out[nod_out_n_col*peconn[ir]+ic] = val;
      }
    }
    peconn += econn_n_col;
    pconn += conn_n_col;
  }

  return( RET_OK );
}

/*!
  For 3D and 2D.

  @par Revision history:
  - 08.06.2006, c
  - 02.08.2006
*/
int32 orient_elements( int32 *flag, int32 flag_n_row,
		      int32 *conn, int32 conn_n_row, int32 conn_n_col,
		      float64 *coors, int32 coors_n_row, int32 coors_n_col,
		      int32 *v_roots, int32 v_roots_n_row,
		      int32 *v_vecs, int32 v_vecs_n_row, int32 v_vecs_n_col,
		      int32 *swap_from, int32 swap_from_n_row, int32 swap_from_n_col,
		      int32 *swap_to, int32 swap_to_n_row, int32 swap_to_n_col )
{
#define IR( iel, ir ) (conn[conn_n_col*(iel)+v_roots[ir]])
#define IV( iel, ir, iv ) (conn[conn_n_col*(iel)+v_vecs[v_vecs_n_col*ir+iv]])
#define CONN( iel, ip ) (conn[conn_n_col*(iel)+ip])
#define SWF( ir, is ) (swap_from[swap_from_n_col*ir+is])
#define SWT( ir, is ) (swap_to[swap_to_n_col*ir+is])

  int32 ir, iel, ii, ip0, ip1, ip2, ip3, tmp, nc;
  float64 v0[3], v1[3], v2[3], v3[3], cross[3], dot[1];

  nc = coors_n_col;
  if (nc == 4) { // 3D.
    for (iel = 0; iel < conn_n_row; iel++) {
      flag[iel] = 0;

      for (ir = 0; ir < v_roots_n_row; ir++) {
	ip0 = IR( iel, ir );
	ip1 = IV( iel, ir, 0 );
	ip2 = IV( iel, ir, 1 );
	ip3 = IV( iel, ir, 2 );
	for (ii = 0; ii < 3; ii++) {
	  v0[ii] = coors[nc*ip0+ii];
	  v1[ii] = coors[nc*ip1+ii] - v0[ii];
	  v2[ii] = coors[nc*ip2+ii] - v0[ii];
	  v3[ii] = coors[nc*ip3+ii] - v0[ii];
	}
	gtr_cross_product( cross, v1, v2 );
	gtr_dot_v3( dot, v3, cross );
/*       output( "%d %d -> %d %d %d %d %e\n", iel, ir, ip0, ip1, ip2, ip3, */
/* 	      dot[0] ); */
	if (dot[0] < CONST_MachEps) {
	  flag[iel]++;
	  for (ii = 0; ii < swap_from_n_col; ii++) {
	    SwapValues( CONN( iel, SWF( ir, ii ) ),
			CONN( iel, SWT( ir, ii ) ), tmp );
/* 	  output( "%d %d\n", SWF( ir, ii ), SWT( ir, ii ) ); */
	  }
	}
      }
/*     sys_pause(); */
    }
  } else if (nc == 3) { // 2D.
    for (iel = 0; iel < conn_n_row; iel++) {
      flag[iel] = 0;

      for (ir = 0; ir < v_roots_n_row; ir++) {
	ip0 = IR( iel, ir );
	ip1 = IV( iel, ir, 0 );
	ip2 = IV( iel, ir, 1 );
	for (ii = 0; ii < 2; ii++) {
	  v0[ii] = coors[nc*ip0+ii];
	  v1[ii] = coors[nc*ip1+ii] - v0[ii];
	  v2[ii] = coors[nc*ip2+ii] - v0[ii];
	}
	v1[2] = v2[2] = 0.0;
	gtr_cross_product( cross, v1, v2 );
	if (cross[2] < CONST_MachEps) {
	  flag[iel]++;
	  for (ii = 0; ii < swap_from_n_col; ii++) {
	    SwapValues( CONN( iel, SWF( ir, ii ) ),
			CONN( iel, SWT( ir, ii ) ), tmp );
	  }
	}
      }
    }
  }
  
  return( RET_OK );

#undef IR
#undef IV
#undef CONN
#undef SWF
#undef SWT
}

/*!
  @par Revision history:
  - 06.03.2005, c
  - 07.03.2005
  - 29.08.2007, from gr_components() (00.01.18)
*/
int32 graph_components( int32 *p_n_comp,
		       int32 *flag, int32 flag_len,
		       int32 *row, int32 row_len,
		       int32 *col, int32 col_len,
		       int32 *pos, int32 pos_len )
{
  // pos is a work array: list of nodes (rows) to process.
  int32 ret = RET_OK, n_tot, n_pos, n_pos_new, n_pos0, n_new, n_stop, n_nod;
  int32 icomp, ii, ir, ic;

  n_nod = row_len - 1;

  n_stop = n_nod;
  for (ir = 0; ir < n_nod; ir++) {
    flag[ir] = -1;
    if ((row[ir+1] - row[ir]) == 0) n_stop--;
  }

  n_tot = 0;
  for (icomp = 0; icomp < n_nod; icomp++) {
    // Find seed.
    ii = 0;
    while (flag[ii] >= 0) {
      ii++;
      if (ii >= n_nod) {
	errput( "error in graph_components()!\n" );
	ERR_CheckGo( ret );
      }
    }
    flag[ii] = icomp;
    pos[0] = ii;
    n_pos0 = 0;
    n_pos_new = n_pos = 1;

    for (ii = 0; ii < n_nod; ii++) {
      n_new = 0;
      for (ir = n_pos0; ir < n_pos; ir++) {
	for (ic = row[pos[ir]]; ic < row[pos[ir]+1]; ic++) {
	  if (flag[col[ic]] == -1) {
	    flag[col[ic]] = icomp;
	    pos[n_pos_new] = col[ic];
	    n_pos_new++;
	    n_new++;
	  }
	}
      }
/*       output( "%d %d %d\n", ii, n_pos, n_new ); */
      n_pos0 = n_pos;
      n_pos = n_pos_new;
      if (n_new == 0) break;
    }
    n_tot += n_pos;
/*     output( "  %d %d %d %d\n", icomp, n_tot, n_stop, n_nod ); */
/*     getchar(); */
    if (n_tot == n_stop) {
      *p_n_comp = icomp + 1;
      break;
    }
  }

 end_label:
  return( ret );
}

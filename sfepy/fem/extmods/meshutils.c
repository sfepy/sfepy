#include <stdlib.h>

#include "meshutils.h"
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

#define SwapValues( a, b, work ) do {\
  (work) = (a); (a) = (b); (b) = (work);\
} while (0)

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
  if (nc == 3) { // 3D.
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
  } else if (nc == 2) { // 2D.
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

#include <stdlib.h>

#include "meshutils.h"
#include "sort.h"
#include "geomtrans.h"

int32 *ic;
int32 nC;

int compareRowCols( const void *objA, const void *objB )
{
  int32 *rowA, *rowB;
  int32 ii;
  int ret = 0;

  rowA = (int32 *) objA;
  rowB = (int32 *) objB;
  
/*   output( "%p %d\n", ic, nC ); */
  for (ii = 0; ii < nC; ii++) {
/*     output( "%d %d\n", rowA[ic[ii]], rowB[ic[ii]]  ); */
/*     output( "%d\n", ic[ii] ); */
    if (rowA[ic[ii]] < rowB[ic[ii]]) {
      ret = -1;
      break;
    } else if (rowA[ic[ii]] > rowB[ic[ii]]) {
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
int32 sortRows( int32 *array, int32 nRow, int32 nCol,
		int32 *iSortCol, int32 nSortCol )
{
/*   output( "%p %p\n", array, iSortCol ); */
  if ((array == 0) || (iSortCol == 0)) {
    output( "null pointer: %p %p\n", array, iSortCol );
  }

  int32_sortRows( array, nRow, nCol, iSortCol, nSortCol );
/*   ic = iSortCol; */
/*   nC = nSortCol; */

/*   qsort( array, nRow, nCol * sizeof( int32 ), &compareRowCols ); */

/*   { */
/*     int32 ir, ic; */

/*     for (ir = 0; ir < nRow; ir++) { */
/*       for (ic = 0; ic < nCol; ic++) { */
/* 	output( "%d ", array[nCol*ir+ic] ); */
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
int32 createList( int32 *p_ii, int32 *objs, int32 objs_nRow, int32 objs_nCol,
		  int32 ig,
		  int32 *conn0, int32 conn_nRow, int32 conn_nCol,
		  int32 *items, int32 items_nRow, int32 items_nCol,
		  int32 isSort )
{
  int32 *conn;
  int32 ii0, ii, ie, ij, io, iop, aux;

/*   output( "%d: %d %d %d %d %d %d\n", *p_ii, objs_nRow, objs_nCol, */
/* 	  conn_nRow, conn_nCol, items_nRow, items_nCol ); */
/*   output( "%p\n", conn0 ); */

  conn = conn0;
  ii0 = ii = *p_ii;
/*   output( "%d\n", ii ); */
  for (ie = 0; ie < conn_nRow; ie++) {
    for (ij = ii; ij < (ii + items_nRow); ij++) {
      io = ij - ii;
      objs[objs_nCol * ij + 0] = ig;
      objs[objs_nCol * ij + 1] = ie;
      objs[objs_nCol * ij + 2] = io;

      for (iop = 0; iop < items_nCol; iop++) {
/* 	output( "%d %d %d\n", ij, iop, conn[items[items_nCol*io+iop]] ); */
	aux = items[items_nCol*io+iop];
	if (aux >= 0)
	  objs[objs_nCol * ij + 3 + iop] = conn[aux];
      }
    }
    ii += items_nRow;
    conn += conn_nCol;
  }
  *p_ii = ii;
/*   output( "%d\n", ii ); */

  if (isSort) {
    objs += ii0 * objs_nCol;
    if (items_nCol == 2) {
      for (ii = ii0; ii < (*p_ii); ii++) {
	MESH_Sort2( objs + 3, aux );
	objs += objs_nCol;
      }
    } else if (items_nCol == 3) {
      for (ii = ii0; ii < (*p_ii); ii++) {
	MESH_Sort3( objs + 3, aux );
	objs += objs_nCol;
      }
    } else if (items_nCol == 4) {
      for (ii = ii0; ii < (*p_ii); ii++) {
	if (objs[6] == -1) {
	  MESH_Sort3( objs + 3, aux );
	} else {
	  MESH_Sort4( objs + 3, aux );
	}
	objs += objs_nCol;
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
int32 neighbourListPtr( int32 *pobj, int32 pobj_nRow,
			int32 *pg, int32 pg_nRow,
			int32 *pel, int32 pel_nRow,
			int32 *dataS, int32 dataS_nRow, int32 dataS_nCol,
			int32 *ic, int32 ic_nRow,
			int32 mode )
{
  int32 ii, i1, ip, icmp, cmp, nCmp;
  int32 *aux0, *aux1;

  if (pobj_nRow == 0) return( RET_OK );

  nCmp = (mode == 0) ? 2 : 4;

  pobj[0] = 0;
  ii = 0;
  while (ii < (dataS_nRow - 2)) {
    if (ic[ii]) {
      ip = ii;
      aux0 = dataS + dataS_nCol * ii + 3;
      while (1) {
	aux1 = dataS + dataS_nCol * (ip + 1) + 3;
	cmp = 1;
	for (icmp = 0; icmp < nCmp; icmp++) {
	  if (aux0[icmp] != aux1[icmp]) {
	    cmp = 0;
	    break;
	  }
	}
	if (!cmp) break;
	ip += 1;
      }
/*       output( "neighbourListPtr %d %d\n", ii, ip ); */
      for (i1 = ii; i1 < (ip + 1); i1++) {
	aux1 = dataS + dataS_nCol * i1;
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
int32 neighbourList( int32 *p_iu, int32 *objs, int32 objs_nRow,
		     int32 *uid, int32 uid_nRow,
		     int32 *cnt, int32 cnt_nRow,
		     int32 *pg, int32 pg_nRow,
		     int32 *pel, int32 pel_nRow,
		     int32 *pobj, int32 pobj_nRow,
		     int32 *dataS, int32 dataS_nRow, int32 dataS_nCol,
		     int32 *uidIn, int32 uidIn_nRow,
		     int32 *ic, int32 ic_nRow,
		     int32 *perm, int32 perm_nRow,
		     int32 mode )
{
  int32 ii, iu, i1, i2, ip, icmp, cmp, nCmp, ptr;
  int32 *aux0, *aux1;

  if (objs_nRow == 0) return( RET_OK );

  nCmp = (mode == 0) ? 2 : 4;

  ii = 0;
  iu = 0;
  while (ii < (dataS_nRow - 2)) {
    if (ic[ii]) {
      ip = ii;
      aux0 = dataS + dataS_nCol * ii + 3;
      while (1) {
	aux1 = dataS + dataS_nCol * (ip + 1) + 3;
	cmp = 1;
	for (icmp = 0; icmp < nCmp; icmp++) {
	  if (aux0[icmp] != aux1[icmp]) {
	    cmp = 0;
	    break;
	  }
	}
	if (!cmp) break;
	ip += 1;
      }
/*       output( "neighbourList %d %d\n", ii, ip ); */
      for (i1 = ii; i1 < (ip + 1); i1++) {
	// Own edge/face first.
	aux1 = dataS + dataS_nCol * i1;
	ptr = pel[pg[aux1[0]] + aux1[1]] + aux1[2];
	objs[pobj[ptr]+cnt[ptr]] = perm[i1];
	cnt[ptr] += 1;
	uid[ptr] = iu;
	if (iu != uidIn[i1]) {
	  errput( "%d == %d\n", iu, uidIn[i1] );
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
      aux1 = dataS + dataS_nCol * ii;
      ptr = pel[pg[aux1[0]] + aux1[1]] + aux1[2];
      objs[pobj[ptr]+cnt[ptr]] = perm[ii];
      cnt[ptr] += 1;
      uid[ptr] = iu;
      if (iu != uidIn[ii]) {
	errput( "%d == %d\n", iu, uidIn[ii] );
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
int32 orientEdges( int32 *ori, int32 ori_nRow, int32 ori_nCol,
		   int32 *conn, int32 conn_nRow, int32 conn_nCol,
		   int32 *edges, int32 edges_nRow, int32 edges_nCol )
{
  int32 ir, ic;

  for (ir = 0; ir < ori_nRow; ir++) {
    for (ic = 0; ic < ori_nCol; ic++) {
      ori[ic] = conn[edges[2*ic+1]] < conn[edges[2*ic]];
    }
    ori += ori_nCol;
    conn += conn_nCol;
  }

  return( RET_OK );
}

/*!
  No shape checks!

  @par Revision history:
  - 10.10.2005, c
*/
int32 assignEdgeNodes( int32 *p_iseq,  
		       int32 *econn, int32 econn_nRow, int32 econn_nCol,
		       int32 *cntEN, int32 cntEN_nRow, int32 cntEN_nCol,
		       int32 *ori, int32 ori_nRow, int32 ori_nCol,
		       int32 *ntt, int32 ntt_nRow, int32 ntt_nCol,
		       int32 *uid, int32 uid_nRow,
		       int32 **edgeDesc, int32 *edgeDesc_nRow, int32 nEdge,
		       int32 cptr0 )
{
  int32 iel, ii, en, iep, it, ui, iseq, nEN;
  int32 *ei, *pori, *puid, *peconn, *pcnt;

  iseq = *p_iseq;
  pori = ori;
  puid = uid + cptr0;
  peconn = econn;
/*   output( "%d %d %d %d %d %d %d %d %d %d %d %d\n", */
/* 	  econn_nRow, econn_nCol, cntEN_nRow, cntEN_nCol, */
/* 	  ori_nRow, ori_nCol, ntt_nRow, ntt_nCol, */
/* 	  uid_nRow, nEdge, cptr0, iseq ); */
  for (iel = 0; iel < econn_nRow; iel++) {

    for (ii = 0; ii < nEdge; ii++) {
      ei = edgeDesc[ii];
      nEN = edgeDesc_nRow[ii];
      for (en = 0; en < nEN; en++) {
	iep = ei[2*en];
	it = ntt[ntt_nCol*pori[ii]+en];
	ui = puid[ii];
	pcnt = cntEN + cntEN_nCol * it + ui;
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
    pori += ori_nCol;
    puid += nEdge;
    peconn += econn_nCol;
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
int32 interpVertexData( float64 *nodOut, int32 nodOut_nRow, int32 nodOut_nCol,
			int32 *econn, int32 econn_nRow, int32 econn_nCol,
			float64 *nodIn, int32 nodIn_nRow, int32 nodIn_nCol,
			int32 *conn, int32 conn_nRow, int32 conn_nCol,
			float64 *bf, int32 bf_nRow, int32 bf_nCol,
			int32 omitCols )
{
  int32 iel, ir, ic, ii;
  int32 *pconn = conn, *peconn = econn;
  float64 val;

  for (iel = 0; iel < econn_nRow; iel++) {
    for (ir = 0; ir < bf_nRow; ir++){
      for (ic = 0; ic < (nodOut_nCol - omitCols); ic++) {
	val = 0.0;
	for (ii = 0; ii < bf_nCol; ii++) {
	  val += bf[bf_nCol*ir+ii] * nodIn[nodIn_nCol*pconn[ii]+ic];
	}
	nodOut[nodOut_nCol*peconn[ir]+ic] = val;
      }
    }
    peconn += econn_nCol;
    pconn += conn_nCol;
  }

  return( RET_OK );
}

/*!
  For 3D and 2D.

  @par Revision history:
  - 08.06.2006, c
  - 02.08.2006
*/
int32 orientElements( int32 *flag, int32 flag_nRow,
		      int32 *conn, int32 conn_nRow, int32 conn_nCol,
		      float64 *coors, int32 coors_nRow, int32 coors_nCol,
		      int32 *vRoots, int32 vRoots_nRow,
		      int32 *vVecs, int32 vVecs_nRow, int32 vVecs_nCol,
		      int32 *swapFrom, int32 swapFrom_nRow, int32 swapFrom_nCol,
		      int32 *swapTo, int32 swapTo_nRow, int32 swapTo_nCol )
{
#define IR( iel, ir ) (conn[conn_nCol*(iel)+vRoots[ir]])
#define IV( iel, ir, iv ) (conn[conn_nCol*(iel)+vVecs[vVecs_nCol*ir+iv]])
#define CONN( iel, ip ) (conn[conn_nCol*(iel)+ip])
#define SWF( ir, is ) (swapFrom[swapFrom_nCol*ir+is])
#define SWT( ir, is ) (swapTo[swapTo_nCol*ir+is])

  int32 ir, iel, ii, ip0, ip1, ip2, ip3, tmp, nc;
  float64 v0[3], v1[3], v2[3], v3[3], cross[3], dot[1];

  nc = coors_nCol;
  if (nc == 4) { // 3D.
    for (iel = 0; iel < conn_nRow; iel++) {
      flag[iel] = 0;

      for (ir = 0; ir < vRoots_nRow; ir++) {
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
	gtr_crossProduct( cross, v1, v2 );
	gtr_dotV3( dot, v3, cross );
/*       output( "%d %d -> %d %d %d %d %e\n", iel, ir, ip0, ip1, ip2, ip3, */
/* 	      dot[0] ); */
	if (dot[0] < CONST_MachEps) {
	  flag[iel]++;
	  for (ii = 0; ii < swapFrom_nCol; ii++) {
	    SwapValues( CONN( iel, SWF( ir, ii ) ),
			CONN( iel, SWT( ir, ii ) ), tmp );
/* 	  output( "%d %d\n", SWF( ir, ii ), SWT( ir, ii ) ); */
	  }
	}
      }
/*     sys_pause(); */
    }
  } else if (nc == 3) { // 2D.
    for (iel = 0; iel < conn_nRow; iel++) {
      flag[iel] = 0;

      for (ir = 0; ir < vRoots_nRow; ir++) {
	ip0 = IR( iel, ir );
	ip1 = IV( iel, ir, 0 );
	ip2 = IV( iel, ir, 1 );
	for (ii = 0; ii < 2; ii++) {
	  v0[ii] = coors[nc*ip0+ii];
	  v1[ii] = coors[nc*ip1+ii] - v0[ii];
	  v2[ii] = coors[nc*ip2+ii] - v0[ii];
	}
	v1[2] = v2[2] = 0.0;
	gtr_crossProduct( cross, v1, v2 );
	if (cross[2] < CONST_MachEps) {
	  flag[iel]++;
	  for (ii = 0; ii < swapFrom_nCol; ii++) {
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
int32 graphComponents( int32 *p_nComp,
		       int32 *flag, int32 flag_len,
		       int32 *row, int32 row_len,
		       int32 *col, int32 col_len,
		       int32 *pos, int32 pos_len )
{
  // pos is a work array: list of nodes (rows) to process.
  int32 ret = RET_OK, nTot, nPos, nPosNew, nPos0, nNew, nStop, nNod;
  int32 icomp, ii, ir, ic;

  nNod = row_len - 1;

  nStop = nNod;
  for (ir = 0; ir < nNod; ir++) {
    flag[ir] = -1;
    if ((row[ir+1] - row[ir]) == 0) nStop--;
  }

  nTot = 0;
  for (icomp = 0; icomp < nNod; icomp++) {
    // Find seed.
    ii = 0;
    while (flag[ii] >= 0) {
      ii++;
      if (ii >= nNod) {
	errput( "error in graphComponents()!\n" );
	ERR_CheckGo( ret );
      }
    }
    flag[ii] = icomp;
    pos[0] = ii;
    nPos0 = 0;
    nPosNew = nPos = 1;

    for (ii = 0; ii < nNod; ii++) {
      nNew = 0;
      for (ir = nPos0; ir < nPos; ir++) {
	for (ic = row[pos[ir]]; ic < row[pos[ir]+1]; ic++) {
	  if (flag[col[ic]] == -1) {
	    flag[col[ic]] = icomp;
	    pos[nPosNew] = col[ic];
	    nPosNew++;
	    nNew++;
	  }
	}
      }
/*       output( "%d %d %d\n", ii, nPos, nNew ); */
      nPos0 = nPos;
      nPos = nPosNew;
      if (nNew == 0) break;
    }
    nTot += nPos;
/*     output( "  %d %d %d %d\n", icomp, nTot, nStop, nNod ); */
/*     getchar(); */
    if (nTot == nStop) {
      *p_nComp = icomp + 1;
      break;
    }
  }

 end_label:
  return( ret );
}

/* -*- C -*- */
#ifdef SWIGPYTHON

%module meshutils

%{
#include "meshutils.h"
%}

%include "types.h"

%include common.i
%include array.i

%apply (int32 *array, int32 nRow, int32 nCol) {
    (int32 *array, int32 nRow, int32 nCol)
};

%apply (int32 *array, int32 len) {
    (int32 *iSortCol, int32 nSortCol)
};

int32 sortRows( int32 *array, int32 nRow, int32 nCol,
		int32 *iSortCol, int32 nSortCol );


%apply (int32 *p_int32_in) {
    (int32 *p_ii)
};
%apply (int32 *p_int32_out) {
    (int32 *p_ii)
};
%apply (int32 *array, int32 nRow, int32 nCol) {
    (int32 *objs, int32 objs_nRow, int32 objs_nCol),
    (int32 *conn0, int32 conn_nRow, int32 conn_nCol),
    (int32 *items, int32 items_nRow, int32 items_nCol)
};
int32 createList( int32 *p_ii, int32 *objs, int32 objs_nRow, int32 objs_nCol,
		  int32 ig,
		  int32 *conn0, int32 conn_nRow, int32 conn_nCol,
		  int32 *items, int32 items_nRow, int32 items_nCol,
		  int32 isSort );

%apply (int32 *array, int32 len) {
    (int32 *pobj, int32 pobj_nRow),
    (int32 *pg, int32 pg_nRow),
    (int32 *pel, int32 pel_nRow),
    (int32 *ic, int32 ic_nRow)
};
%apply (int32 *array, int32 nRow, int32 nCol) {
    (int32 *dataS, int32 dataS_nRow, int32 dataS_nCol)
};
int32 neighbourListPtr( int32 *pobj, int32 pobj_nRow,
			int32 *pg, int32 pg_nRow,
			int32 *pel, int32 pel_nRow,
			int32 *dataS, int32 dataS_nRow, int32 dataS_nCol,
			int32 *ic, int32 ic_nRow,
			int32 mode );


%apply (int32 *p_int32_argout) {
    (int32 *p_iu)
};
%apply (int32 *array, int32 len) {
    (int32 *objs, int32 objs_nRow),
    (int32 *uid, int32 uid_nRow),
    (int32 *uidIn, int32 uidIn_nRow),
    (int32 *cnt, int32 cnt_nRow),
    (int32 *perm, int32 perm_nRow)
};
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
		     int32 mode );

%apply (int32 *array, int32 nRow, int32 nCol) {
    (int32 *ori, int32 ori_nRow, int32 ori_nCol),
    (int32 *conn, int32 conn_nRow, int32 conn_nCol),
    (int32 *edges, int32 edges_nRow, int32 edges_nCol)
};
int32 orientEdges( int32 *ori, int32 ori_nRow, int32 ori_nCol,
		   int32 *conn, int32 conn_nRow, int32 conn_nCol,
		   int32 *edges, int32 edges_nRow, int32 edges_nCol );


%apply (int32 *p_int32_in) {
    (int32 *p_iseq)
};
%apply (int32 *p_int32_out) {
    (int32 *p_iseq)
};

%apply (int32 *array, int32 nRow, int32 nCol) {
    (int32 *econn, int32 econn_nRow, int32 econn_nCol),
    (int32 *cntEN, int32 cntEN_nRow, int32 cntEN_nCol),
    (int32 *ori, int32 ori_nRow, int32 ori_nCol),
    (int32 *ntt, int32 ntt_nRow, int32 ntt_nCol)
};

%apply (int32 *array, int32 len) {
    (int32 *uid, int32 uid_nRow)
};

/*!
  @par Revision history:
  - 10.10.2005, c
*/
%typemap( in ) (int32 **itemList, int32 *item_nRow, int32 nItem) {
  if (PyList_Check( $input )) {
    int32 ii;
    PyArrayObject *obj;
    $3 = PyList_Size( $input );
    
    $1 = allocMem( int32 *, $3 );
    $2 = allocMem( int32, $3 );
    for (ii = 0; ii < $3; ii++) {
      PyObject *oo = PyList_GetItem( $input, ii );
      obj = helper_getCArrayObject( oo, PyArray_INT32, 2, 2 );
      if (!obj) return NULL;
      $1[ii] = (int32 *) obj->data;
      $2[ii] = obj->dimensions[0];
    }
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}
%typemap( freearg ) (int32 **itemList, int32 *item_nRow, int32 nItem) {
  freeMem( (int32 **) $1 );
  freeMem( (int32 *) $2 );
}

%apply (int32 **itemList, int32 *item_nRow, int32 nItem) {
  (int32 **edgeDesc, int32 *edgeDesc_nRow, int32 nEdge)
};

int32 assignEdgeNodes( int32 *p_iseq,  
		       int32 *econn, int32 econn_nRow, int32 econn_nCol,
		       int32 *cntEN, int32 cntEN_nRow, int32 cntEN_nCol,
		       int32 *ori, int32 ori_nRow, int32 ori_nCol,
		       int32 *ntt, int32 ntt_nRow, int32 ntt_nCol,
		       int32 *uid, int32 uid_nRow,
		       int32 **edgeDesc, int32 *edgeDesc_nRow, int32 nEdge,
		       int32 cptr0 );

%apply (int32 *array, int32 nRow, int32 nCol) {
    (int32 *econn, int32 econn_nRow, int32 econn_nCol),
    (int32 *conn, int32 conn_nRow, int32 conn_nCol)
};

%apply (float64 *array, int32 nRow, int32 nCol) {
    (float64 *nodOut, int32 nodOut_nRow, int32 nodOut_nCol),
    (float64 *nodIn, int32 nodIn_nRow, int32 nodIn_nCol),
    (float64 *bf, int32 bf_nRow, int32 bf_nCol)
};

int32 interpVertexData( float64 *nodOut, int32 nodOut_nRow, int32 nodOut_nCol,
			int32 *econn, int32 econn_nRow, int32 econn_nCol,
			float64 *nodIn, int32 nodIn_nRow, int32 nodIn_nCol,
			int32 *conn, int32 conn_nRow, int32 conn_nCol,
			float64 *bf, int32 bf_nRow, int32 bf_nCol,
			int32 omitCols );

%apply (float64 *array, int32 nRow, int32 nCol) {
    (float64 *coors, int32 coors_nRow, int32 coors_nCol)
};

%apply (int32 *array, int32 len) {
    (int32 *flag, int32 flag_nRow),
    (int32 *vRoots, int32 vRoots_nRow)
};

%apply (int32 *array, int32 nRow, int32 nCol) {
    (int32 *conn, int32 conn_nRow, int32 conn_nCol),
    (int32 *vVecs, int32 vVecs_nRow, int32 vVecs_nCol),
    (int32 *swapFrom, int32 swapFrom_nRow, int32 swapFrom_nCol),
    (int32 *swapTo, int32 swapTo_nRow, int32 swapTo_nCol)
};

int32 orientElements( int32 *flag, int32 flag_nRow,
		      int32 *conn, int32 conn_nRow, int32 conn_nCol,
		      float64 *coors, int32 coors_nRow, int32 coors_nCol,
		      int32 *vRoots, int32 vRoots_nRow,
		      int32 *vVecs, int32 vVecs_nRow, int32 vVecs_nCol,
		      int32 *swapFrom, int32 swapFrom_nRow, int32 swapFrom_nCol,
		      int32 *swapTo, int32 swapTo_nRow, int32 swapTo_nCol );

%apply (int32 *array, int32 len) {
    (int32 *flag, int32 flag_len),
    (int32 *row, int32 row_len),
    (int32 *col, int32 col_len),
    (int32 *pos, int32 pos_len)
};
%apply (int32 *p_int32_argout) {
    (int32 *p_nComp)
};
int32 graphComponents( int32 *p_nComp,
		       int32 *flag, int32 flag_len,
		       int32 *row, int32 row_len,
		       int32 *col, int32 col_len,
		       int32 *pos, int32 pos_len );

// %{
// PyObject *isSequence( PyObject *input ) {

//   if (PySequence_Check( input )) {
//     return( PyBool_FromLong( 1 ) );
//   } else {
//     return( PyBool_FromLong( 0 ) );
//   }
// }
// %}

// PyObject *isSequence( PyObject *input );

#endif

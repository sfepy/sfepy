/* -*- C -*- */
#ifdef SWIGPYTHON

%module fem

%{
#include "fem.h"
%}

%include "types.h"

%include common.i
%include array.i
%include fmfield.i

%apply (FMField *in) {
  (FMField *vec),
  (FMField *vecInEls),
  (FMField *vec_r),
  (FMField *vecInEls_r),
  (FMField *vec_i),
  (FMField *vecInEls_i),
  (FMField *mtx),
  (FMField *mtxInEls),
  (FMField *mtx_r),
  (FMField *mtxInEls_r),
  (FMField *mtx_i),
  (FMField *mtxInEls_i)
};
%apply (int32 *array, int32 n_row, int32 n_col) {
  (int32 *conn, int32 nEl, int32 nEP),
  (int32 *connR, int32 nElR, int32 nEPR),
  (int32 *connC, int32 nElC, int32 nEPC)
};
%apply (int32 *array, int32 len) {
  (int32 *iels, int32 iels_len),
  (int32 *eq, int32 nEq),
  (int32 *prows, int32 prows_len),
  (int32 *cols, int32 cols_len)
};

int32 assemble_vector( FMField *vec, FMField *vecInEls,
		       int32 *iels, int32 iels_len,
		       float64 sign, int32 *conn, int32 nEl, int32 nEP );

int32 assemble_vector_complex( FMField *vec_r, FMField *vec_i,
			       FMField *vecInEls_r, FMField *vecInEls_i,
			       int32 *iels, int32 iels_len,
			       float64 sign_r, float64 sign_i,
			       int32 *conn, int32 nEl, int32 nEP );

int32 assemble_matrix( FMField *mtx,
		      int32 *prows, int32 prows_len,
		      int32 *cols, int32 cols_len,
		      FMField *mtxInEls,
		      int32 *iels, int32 iels_len, float64 sign,
		      int32 *connR, int32 nElR, int32 nEPR,
		      int32 *connC, int32 nElC, int32 nEPC );

int32 assemble_matrix_complex( FMField *mtx_r, FMField *mtx_i,
			       int32 *prows, int32 prows_len,
			       int32 *cols, int32 cols_len,
			       FMField *mtxInEls_r, FMField *mtxInEls_i,
			       int32 *iels, int32 iels_len,
			       float64 sign_r, float64 sign_i,
			       int32 *connR, int32 nElR, int32 nEPR,
			       int32 *connC, int32 nElC, int32 nEPC );

/*!
  @par Revision history:
  - 03.03.2005, c
*/
%typemap( in ) (int32 *nEl, int32 *nEP, int32 **conn) {
  PyObject *aux;
  PyArrayObject *obj;
  int32 ig, nGr;
  int32 *tnEP, *tnEl;
  int32 **tconn;

  if (!PyList_Check( $input )) {
    PyErr_SetString( PyExc_TypeError, "not a list" );
    return NULL;
  }

  nGr = PyList_Size( $input );
  tnEl = alloc_mem( int32, nGr );
  tnEP = alloc_mem( int32, nGr );
  tconn = alloc_mem( int32 *, nGr );
  for (ig = 0; ig < nGr; ig++) {
    aux = PyList_GetItem( $input, ig );
    obj = helper_get_c_array_object( aux, PyArray_INT32, 0, 0 );
    if (!obj) return NULL;
    tnEl[ig] = obj->dimensions[0];
    tnEP[ig] = obj->dimensions[1];
    tconn[ig] = (int32 *) obj->data;
    Py_DECREF( obj );
  }

  $1 = tnEl;
  $2 = tnEP;
  $3 = tconn;
};
%typemap( freearg ) (int32 *nEl, int32 *nEP, int32 **conn) {
  free_mem( $1 );
  free_mem( $2 );
  free_mem( $3 );
}

%apply (int32 *nEl, int32 *nEP, int32 **conn) {
  (int32 *nElR, int32 *nEPR, int32 **connR),
  (int32 *nElC, int32 *nEPC, int32 **connC)
};

%apply (int32 *p_len, int32 **p_array) {
  (int32 *p_nRow, int32 **p_prow),
  (int32 *p_nnz, int32 **p_icol)
};

int32 raw_graph( int32 *p_nRow, int32 **p_prow,
		 int32 *p_nnz, int32 **p_icol,
		 int32 nRow, int32 nCol, int32 nGr,
		 int32 *nElR, int32 *nEPR, int32 **connR,
		 int32 *nElC, int32 *nEPC, int32 **connC );

#endif

/* -*- C -*- */
#ifdef SWIGPYTHON

%module meshutils

%{
#include "meshutils.h"
%}

%include "types.h"

%include common.i
%include array.i

%apply (int32 *p_int32_in) {
    (int32 *p_iseq)
};
%apply (int32 *p_int32_out) {
    (int32 *p_iseq)
};

%apply (int32 *array, int32 n_row, int32 n_col) {
    (int32 *econn, int32 econn_n_row, int32 econn_n_col),
    (int32 *cnt_en, int32 cnt_en_n_row, int32 cnt_en_n_col),
    (int32 *ori, int32 ori_n_row, int32 ori_n_col),
    (int32 *ntt, int32 ntt_n_row, int32 ntt_n_col)
};

%apply (int32 *array, int32 len) {
    (int32 *uid, int32 uid_n_row)
};

/*!
  @par Revision history:
  - 10.10.2005, c
*/
%typemap( in ) (int32 **item_list, int32 *item_n_row, int32 n_item) {
  if (PyList_Check( $input )) {
    int32 ii;
    PyArrayObject *obj;
    $3 = PyList_Size( $input );
    
    $1 = alloc_mem( int32 *, $3 );
    $2 = alloc_mem( int32, $3 );
    for (ii = 0; ii < $3; ii++) {
      PyObject *oo = PyList_GetItem( $input, ii );
      obj = helper_get_c_array_object( oo, PyArray_INT32, 1, 1 );
      if (!obj) return NULL;
      $1[ii] = (int32 *) obj->data;
      $2[ii] = obj->dimensions[0];
    }
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}
%typemap( freearg ) (int32 **item_list, int32 *item_n_row, int32 n_item) {
  free_mem( (int32 **) $1 );
  free_mem( (int32 *) $2 );
}

%apply (int32 **item_list, int32 *item_n_row, int32 n_item) {
  (int32 **edge_desc, int32 *edge_desc_n_row, int32 n_edge)
};

int32 assign_edge_nodes( int32 *p_iseq,  
		       int32 *econn, int32 econn_n_row, int32 econn_n_col,
		       int32 *cnt_en, int32 cnt_en_n_row, int32 cnt_en_n_col,
		       int32 *ori, int32 ori_n_row, int32 ori_n_col,
		       int32 *ntt, int32 ntt_n_row, int32 ntt_n_col,
		       int32 *uid, int32 uid_n_row,
		       int32 **edge_desc, int32 *edge_desc_n_row, int32 n_edge,
		       int32 cptr0 );

%apply (int32 *array, int32 n_row, int32 n_col) {
    (int32 *econn, int32 econn_n_row, int32 econn_n_col),
    (int32 *conn, int32 conn_n_row, int32 conn_n_col)
};

%apply (float64 *array, int32 n_row, int32 n_col) {
    (float64 *nod_out, int32 nod_out_n_row, int32 nod_out_n_col),
    (float64 *nod_in, int32 nod_in_n_row, int32 nod_in_n_col),
    (float64 *bf, int32 bf_n_row, int32 bf_n_col)
};

int32 interp_vertex_data( float64 *nod_out, int32 nod_out_n_row, int32 nod_out_n_col,
			int32 *econn, int32 econn_n_row, int32 econn_n_col,
			float64 *nod_in, int32 nod_in_n_row, int32 nod_in_n_col,
			int32 *conn, int32 conn_n_row, int32 conn_n_col,
			float64 *bf, int32 bf_n_row, int32 bf_n_col,
			int32 omit_cols );

%apply (float64 *array, int32 n_row, int32 n_col) {
    (float64 *coors, int32 coors_n_row, int32 coors_n_col)
};

%apply (int32 *array, int32 len) {
    (int32 *flag, int32 flag_n_row),
    (int32 *v_roots, int32 v_roots_n_row)
};

%apply (int32 *array, int32 n_row, int32 n_col) {
    (int32 *conn, int32 conn_n_row, int32 conn_n_col),
    (int32 *v_vecs, int32 v_vecs_n_row, int32 v_vecs_n_col),
    (int32 *swap_from, int32 swap_from_n_row, int32 swap_from_n_col),
    (int32 *swap_to, int32 swap_to_n_row, int32 swap_to_n_col)
};

int32 orient_elements( int32 *flag, int32 flag_n_row,
		      int32 *conn, int32 conn_n_row, int32 conn_n_col,
		      float64 *coors, int32 coors_n_row, int32 coors_n_col,
		      int32 *v_roots, int32 v_roots_n_row,
		      int32 *v_vecs, int32 v_vecs_n_row, int32 v_vecs_n_col,
		      int32 *swap_from, int32 swap_from_n_row, int32 swap_from_n_col,
		      int32 *swap_to, int32 swap_to_n_row, int32 swap_to_n_col );

%apply (int32 *array, int32 len) {
    (int32 *flag, int32 flag_len),
    (int32 *row, int32 row_len),
    (int32 *col, int32 col_len),
    (int32 *pos, int32 pos_len)
};
%apply (int32 *p_int32_argout) {
    (int32 *p_n_comp)
};
int32 graph_components( int32 *p_n_comp,
		       int32 *flag, int32 flag_len,
		       int32 *row, int32 row_len,
		       int32 *col, int32 col_len,
		       int32 *pos, int32 pos_len );

// %{
// PyObject *is_sequence( PyObject *input ) {

//   if (PySequence_Check( input )) {
//     return( PyBool_FromLong( 1 ) );
//   } else {
//     return( PyBool_FromLong( 0 ) );
//   }
// }
// %}

// PyObject *is_sequence( PyObject *input );

#endif

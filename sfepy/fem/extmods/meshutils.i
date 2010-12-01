/* -*- C -*- */
#ifdef SWIGPYTHON

%module meshutils

%{
#include "meshutils.h"
%}

%include "types.h"

%include common.i
%include array.i

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

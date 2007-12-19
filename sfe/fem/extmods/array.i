/* -*- C -*- */
#ifdef SWIGPYTHON

%module array

%{
#include "numpy/arrayobject.h"
  //#include "numpy/noprefix.h"
%}

/*!
  @par Revision history:
  - 20.02.2004, c
*/
%init %{
    import_array();
//     printf( "module 'array' initialized.\n" );
%}

%{
/*!
  @par Revision history:
  - 22.02.2005, c
  - 03.03.2005
  - 25.11.2005
  - 30.11.2005
  - 01.12.2005
*/
PyArrayObject *helper_getCArrayObject( PyObject *input, int type,
				       int minDim, int maxDim ) {
  PyArrayObject *obj;

  if (PyArray_Check( input )) {
    obj = (PyArrayObject *) input;
    if (!PyArray_ISCARRAY( obj )) {
      PyErr_SetString( PyExc_TypeError, "not a C array" );
      return NULL;
    }
    obj = (PyArrayObject *)
      PyArray_ContiguousFromAny( input, type, minDim, maxDim );
    if (!obj) return NULL;
  } else {
    PyErr_SetString( PyExc_TypeError, "not an array" );
    return NULL;
  }
  return obj;
}

/*!
  Makes a copy -> free original buffer afterwards.

  @par Revision history:
  - 24.02.2004, c
  - 27.11.2005, adapted from mafest1
*/
PyArrayObject *helper_newCArrayObject_i32( int32 len, int32 *array ) {
  int32 ii;
  int32 *out;
  npy_intp plen[1];
  PyArrayObject *obj = 0;

  plen[0] = len;
  obj = (PyArrayObject *) PyArray_SimpleNew( 1, plen, PyArray_INT32 );

  if (obj) {
    out = (int32 *) PyArray_DATA( obj );
    for (ii = 0; ii < len; ii++) {
      out[ii] = array[ii];
    }
  } else {
    PyErr_SetString( PyExc_ValueError, "array not created" );
    return NULL;
  }
  return( obj );
}
%}

/*!
  @par Revision history:
  - 03.03.2005, c
*/
%typemap( in ) (int32 *array) {
  PyArrayObject *obj;

  obj = helper_getCArrayObject( $input, PyArray_INT32, 0, 0 );
  if (!obj) return NULL;

  $1 = (int32 *) PyArray_DATA( obj );
  Py_DECREF( obj );
};

/*!
  @par Revision history:
  - 14.12.2004, c
  - 22.02.2005
*/
%typemap( in ) (int32 *array, int32 len) {
  PyArrayObject *obj;

  obj = helper_getCArrayObject( $input, PyArray_INT32, 1, 1 );
  if (!obj) return NULL;

  $1 = (int32 *) PyArray_DATA( obj );
  $2 = PyArray_DIM( obj, 0 );
  Py_DECREF( obj );
};

/*!
  @par Revision history:
  - 14.12.2004, c
  - 22.02.2005
*/
%typemap( in ) (int32 *array, int32 nRow, int32 nCol) {
  PyArrayObject *obj;

  obj = helper_getCArrayObject( $input, PyArray_INT32, 2, 2 );
  if (!obj) return NULL;

  $1 = (int32 *) PyArray_DATA( obj );
  $2 = PyArray_DIM( obj, 0 );
  $3 = PyArray_DIM( obj, 1 );
  Py_DECREF( obj );
};

/*!
  @par Revision history:
  - 10.10.2005, c
*/
%typemap( in ) (float64 *array, int32 nRow, int32 nCol) {
  PyArrayObject *obj;

  obj = helper_getCArrayObject( $input, PyArray_FLOAT64, 2, 2 );
  if (!obj) return NULL;

  $1 = (float64 *) PyArray_DATA( obj );
  $2 = PyArray_DIM( obj, 0 );
  $3 = PyArray_DIM( obj, 1 );
  Py_DECREF( obj );
};

/*!
  @par Revision history:
  - 27.11.2005, c
*/
%typemap( in, numinputs=0 ) (int32 *p_len, int32 **p_array)
     ( int32 tmp1, int32 *tmp2 ) {
  $1 = &tmp1;
  $2 = &tmp2;
}

/*!
  @par Revision history:
  - 27.11.2005, c
*/
%typemap( argout ) (int32 *p_len, int32 **p_array) {
  PyArrayObject *obj;

  obj = helper_newCArrayObject_i32( *$1, *$2 );
  freeMem( *$2 );
  if (obj == NULL) return( NULL );

  $result = helper_appendToTuple( $result, PyArray_Return( obj ) );
}
#endif

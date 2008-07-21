/* -*- C -*- */
#ifdef SWIGPYTHON

%module fmfield

%{
#include "fmfield.h"
%}

%{
/*!
  @par Revision history:
  - 09.12.2005, c
*/
int32 helper_pretend_FMField( FMField *out, PyObject *input )
{
  PyArrayObject *obj;
  int32 ii, nCell = 1, nLev = 1, nRow = 1, nCol = 1;

  obj = helper_get_c_array_object( input, PyArray_FLOAT64, 1, 4 );
  if (!obj) return 0;

  ii = 0;
  switch (obj->nd) {
  case 4:
    nCell = obj->dimensions[ii];
    ii++;
  case 3:
    nLev = obj->dimensions[ii];
    ii++;
  case 2:
    nRow = obj->dimensions[ii];
    ii++;
  case 1:
    nCol = obj->dimensions[ii];
  }
//   printf( "asd %d, %d %d %d %d\n", ii, nCell, nLev, nRow, nCol );

  out->nAlloc = -1;
  fmf_pretend( out, nCell, nLev, nRow, nCol, (float64 *) obj->data );
  Py_DECREF( obj );
  return( 1 );
}
%}

/*!
  @par Revision history:
  - 11.10.2005, c
  - 09.12.2005
*/
%typemap( in ) (FMField *in) (FMField out[1]) {

  if (helper_pretend_FMField( out, $input )) {
    $1 = out;
  } else {
    return NULL;
  }
};

#endif

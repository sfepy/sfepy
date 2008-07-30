/* -*- C -*- */
#ifdef SWIGPYTHON

%module fmfield

%{
#include "fmfield.h"
%}

%{
int32 helper_pretend_FMField( FMField *out, PyObject *input )
{
  PyArrayObject *obj;
  int32 ii, stride, nCell = 1, nLev = 1, nRow = 1, nCol = 1;

  obj = helper_get_c_array_object( input, PyArray_FLOAT64, 1, 4 );
  if (!obj) {
    obj = helper_get_array_object( input, 1, 4 );
    if (!obj) return 0;
    else PyErr_Clear();
  }

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
/*   printf( "dims: %d\n", obj->nd ); */
/*   {int32 jj; */
/*     for (jj = 0; jj < obj->nd; jj++) { */
/*       printf( "%d %d\n", obj->dimensions[jj], obj->strides[jj] ); */
/*     } */
/*   } */
/*   printf( "asd %d, %d %d %d %d\n", ii, nCell, nLev, nRow, nCol ); */

  stride = obj->strides[obj->nd-1];

  out->nAlloc = -1;
  fmf_pretend( out, nCell, nLev, nRow, nCol, (float64 *) obj->data );
  if (stride == 8) { // float64
    // Use offset for stride, fem.i relies on that.
    out->offset = 1;
  } else if (stride == 16) { // real or imag of complex128
    // Use offset for stride, fem.i relies on that.
    out->offset = 2;
    // FMF_SetCell() should work.
    out->cellSize *= 2;
  } else {
    PyErr_SetString( PyExc_TypeError, "unknown array type" );
    return 0;
  }
/*   printf( "%d %d %d\n", out->cellSize, out->offset, out->nColFull ); */

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

/* -*- C -*- */
#ifdef SWIGPYTHON

%module common

%{
#include "common.h"
%}

%feature("autodoc", "1");

/*!
  @par Revision history:
  - 20.02.2004, c
*/
%exception {
  $action
    if (PyErr_Occurred()) {
      return NULL;
    }
}

/*!
  @par Revision history:
  - 17.02.2004, c
*/
%typemap( in ) FILE * {
  if (!PyFile_Check( $input )) {
    PyErr_SetString( PyExc_TypeError, "not a file" );
    return NULL;
  }
  $1 = PyFile_AsFile( $input );
}

/*!
  @par Revision history:
  - 17.02.2005, c
*/
%{
PyObject *helper_append_to_tuple( PyObject *where, PyObject *what ) {
  PyObject *o2, *o3;

  if ((!where) || (where == Py_None)) {
    where = what;
  } else {
    if (!PyTuple_Check( where )) {
      o2 = where;
      where = PyTuple_New( 1 );
      PyTuple_SetItem( where, 0, o2 );
    }
    o3 = PyTuple_New( 1 );
    PyTuple_SetItem( o3, 0, what );
    o2 = where;
    where = PySequence_Concat( o2, o3 );
    Py_DECREF( o2 );
    Py_DECREF( o3 );
  }
  return where;
}
%}

/*!
  @par Revision history:
  - 14.01.2005, c
*/
%typemap( in ) int32 *p_int32_in {
  int32 ii;

  if (PyInt_Check( $input )) {
    ii = PyInt_AsLong( $input );
/*     printf( "p_int32_in: %d\n", ii ); */
  } else {
    PyErr_SetString( PyExc_TypeError, "not an int" );
    return NULL;
  }
  $1 = &ii;
/*   printf( "p_int32_in: %d\n", *$1 ); */
}

/*!
  @par Revision history:
  - 18.01.2005, c
*/
%typemap( in, numinputs=0 ) int32 *p_int32_ignore( int32 tmp ) {
  $1 = &tmp;
}
/*!
  @par Revision history:
  - 14.01.2005, c
  - 17.02.2005
*/
%typemap( argout ) int32 *p_int32_out {
  PyObject *o;

  o = PyInt_FromLong( *$1 );
  $result = helper_append_to_tuple( $result, o );
}

/*!
  @par Revision history:
  - 06.03.2005, c
*/
%typemap( in, numinputs=0 ) int32 *p_int32_argout( int32 tmp ) {
  $1 = &tmp;
}
/*!
  @par Revision history:
  - 06.03.2005, c
*/
%typemap( argout ) int32 *p_int32_argout {
  PyObject *o;

  o = PyInt_FromLong( *$1 );
  $result = helper_append_to_tuple( $result, o );
}

/*!
  @par Revision history:
  - 06.03.2005, c
*/
%typemap( in, numinputs=0 )
     PyObject **p_PyObject_argout( PyObject *tmp ) {
  $1 = &tmp;
}
/*!
  @par Revision history:
  - 06.03.2005, c
*/
%typemap( argout ) PyObject **p_PyObject_argout {
  $result = helper_append_to_tuple( $result, *$1 );
}

void errclear();

void mem_checkIntegrity( int lineNo, char *funName,
			 char *fileName, char *dirName );
void mem_statistics( int lineNo, char *funName,
		     char *fileName, char *dirName );
int32 mem_print( FILE *file, int32 mode );
int32 mem_printSome( FILE *file, int32 mode, int32 num );
int32 mem_freeGarbage();
void sys_pause();

/*!
  @par Revision history:
  - 17.02.2005, c
*/
%inline %{
#define __FUNC__ "??"
void f_checkMemoryIntegrity() {
  mem_checkIntegrity( __LINE__, __FUNC__, __FILE__, __SDIR__ );
}
void f_printMemStats() {
  mem_statistics( __LINE__, __FUNC__, __FILE__, __SDIR__ );
}
%}

#endif

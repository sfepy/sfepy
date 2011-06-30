%module geometry

%{
#include "geometry.h"
%}

%include "types.h"

%include common.i
%include array.i
%include fmfield.i

typedef enum GeometryMode {
  GM_Material,
  GM_Spatial,
} GeometryMode;

typedef struct VolumeGeometry {
  GeometryMode mode;
  int32 nEl;
  int32 nQP;
  int32 dim;
  int32 nEP;
  FMField *bfGM;
  FMField *det; // detJMR or detJSR.
  FMField *volume;
  float64 totalVolume;
} VolumeGeometry;

%extend VolumeGeometry {
  /*!
    @par Revision history:
    - 11.10.2005, c
  */
  VolumeGeometry( int32 nEl, int32 nQP, int32 dim, int32 nEP ) {
    VolumeGeometry *obj;
    /*     printf( "%d %d %d %d\n", nEl, nQP, dim, nEP );*/
    vg_createAlloc( &obj, nEl, nQP, dim, nEP );
    return( obj );
  }

  /*!
    @par Revision history:
    - 11.10.2005, c
  */
  ~VolumeGeometry() {
    vg_freeDestroy( &self );
  }

  /*!
    @par Revision history:
    - 11.10.2005, c
    - 12.10.2006
  */
  char *__str__() {
    static char tmp[256];
    snprintf( tmp, 255,
	      "VolumeGeometry: mode %d, nEl %d, nQP %d, dim: %d, nEP: %d\n",
	      self->mode, self->nEl, self->nQP, self->dim, self->nEP );
    return tmp;
  }

  /*!
    @par Revision history:
    - 11.10.2005, c
  */
  int32 str( FILE *file, int32 mode ) {
    return( vg_print( self, file, mode ) );
  }

  %apply (int32 *array, int32 n_row, int32 n_col) {
      (int32 *conn, int32 nEl, int32 nEP),
      (int32 *edges, int32 edges_nRow, int32 edges_nCol)
  };
  %apply (float64 *array, int32 n_row, int32 n_col) {
      (float64 *coorIn, int32 nNod, int32 dim)
  };
  %apply (FMField *in) {
      (FMField *bfGR),
      (FMField *ebfGR),
      (FMField *weight)
  };

  /*!
    @par Revision history:
    - 11.10.2005, c
  */
  int32 describe( float64 *coorIn, int32 nNod, int32 dim,
		  int32 *conn, int32 nEl, int32 nEP,
		  FMField *bfGR, FMField *ebfGR, FMField *weight ) {
    return( vg_describe( self, coorIn, nNod, dim,
			 conn, nEl, nEP, bfGR, ebfGR, weight ) );
  }

  /*!
    @par Revision history:
    - 12.10.2005, c
    - 25.11.2005
    - 12.04.2007
  */
   PyObject *variable( int32 which ) {
     PyArrayObject *out;
     FMField *obj;
     int32 dims[4];

     switch (which) {
     case 0:
       obj = self->bfGM;
       break;
     case 1:
       obj = self->det;
       break;
     case 2:
       obj = self->volume;
       break;
     default:
       errput( "valid variable range is 0 - 2\n" );
       return 0;
     }
     dims[0] = obj->nCell;
     dims[1] = obj->nLev;
     dims[2] = obj->nRow;
     dims[3] = obj->nCol;
     out = (PyArrayObject *) PyArray_FromDims( 4, dims, PyArray_FLOAT64 );
     memcpy( (float64 *) out->data, obj->val0,
	     obj->nCell * obj->cellSize * sizeof( float64 ) );

     return( (PyObject *) out );
  }

  %apply (FMField *in) {
      (FMField *out),
      (FMField *in)
  };
  %apply (int32 *array, int32 len) {
      (int32 *elList, int32 elList_nRow)
  };

  /*!
    @par Revision history:
    - 15.12.2005, c
  */
   int32 integrate( FMField *out, FMField *in, int32 mode=0 ) {
     return( vg_integrate( self, out, in, mode ) );
  }
  /*!
    @par Revision history:
    - 01.11.2007, c
  */
  int32 integrate_chunk( FMField *out, FMField *in,
			 int32 *elList, int32 elList_nRow ) {
    return( vg_integrateChunk( self, out, in, elList, elList_nRow ) );
  }
  /*!
    @par Revision history:
    - 09.01.2006, c
  */
  int32 get_element_diameters( FMField *out,
			       int32 *edges, int32 edges_nRow, int32 edges_nCol,
			       float64 *coorIn, int32 nNod, int32 dim,
			       int32 *conn, int32 nEl, int32 nEP,
			       int32 *elList, int32 elList_nRow,
			       int32 mode ) {
    return( vg_getElementDiameters( self, out,
				    edges, edges_nRow, edges_nCol,
				    coorIn, nNod, dim,
				    conn, nEl, nEP,
				    elList, elList_nRow,
				    mode ) );
  }
}

typedef struct SurfaceGeometry {
  GeometryMode mode;
  int32 nFa;
  int32 nQP;
  int32 dim;
  int32 nFP;
  FMField *normal;
  FMField *det; // detJMR.
  FMField *bfBGM; // Not computed yet.

  FMField *area;
  float64 totalArea;
} SurfaceGeometry;

%extend SurfaceGeometry {
  /*!
    @par Revision history:
    - 21.12.2005, c
  */
  SurfaceGeometry( int32 nFa, int32 nQP, int32 dim, int32 nFP ) {
    SurfaceGeometry *obj;
    sg_createAlloc( &obj, nFa, nQP, dim, nFP );
    return( obj );
  }

  /*!
    @par Revision history:
    - 21.12.2005, c
  */
  ~SurfaceGeometry() {
    sg_freeDestroy( &self );
  }

  /*!
    @par Revision history:
    - 04.05.2007, c
  */
  int32 alloc_extra_data( int32 nEP  ) {
    fmf_createAlloc( &(self->bfBGM), self->nFa, self->nQP, self->dim, nEP );
    return( RET_OK );
  }

  /*!
    @par Revision history:
    - 21.12.2005, c
  */
  char *__str__() {
    static char tmp[256];
    snprintf( tmp, 255,
	      "SurfaceGeometry: mode %d, nFa %d, nQP %d, dim: %d, nFP: %d",
	      self->mode, self->nFa, self->nQP, self->dim, self->nFP );
    return tmp;
  }

  /*!
    @par Revision history:
    - 21.12.2005, c
  */
  int32 str( FILE *file, int32 mode ) {
    return( sg_print( self, file, mode ) );
  }

  %apply (int32 *array, int32 n_row, int32 n_col) {
      (int32 *fconn, int32 nFa, int32 nFP),
      (int32 *conn, int32 nEl, int32 nEP),
      (int32 *fis, int32 nFa, int32 nFP)
  };
  %apply (float64 *array, int32 n_row, int32 n_col) {
      (float64 *coorIn, int32 nNod, int32 dim)
  };
  %apply (FMField *in) {
      (FMField *bfGR),
      (FMField *bfBGR),
      (FMField *ebfBGR),
      (FMField *weight)
  };

  /*!
    @par Revision history:
    - 21.12.2005, c
  */
  int32 describe( float64 *coorIn, int32 nNod, int32 dim,
		  int32 *fconn, int32 nFa, int32 nFP,
		  FMField *bfGR, FMField *weight ) {
    return( sg_describe( self, coorIn, nNod, dim,
			 fconn, nFa, nFP, bfGR, weight ) );
  }

  /*!
    @par Revision history:
    - 21.12.2005, c
    - 12.04.2007
    - 04.05.2007
  */
  PyObject *variable( int32 which ) {
     PyArrayObject *out;
     FMField *obj;
     int32 dims[4];

     switch (which) {
     case 0:
       obj = self->normal;
       break;
     case 1:
       obj = self->det;
       break;
     case 2:
       obj = self->area;
       break;
     case 3:
       obj = self->bfBGM;
       break;
     default:
       errput( "valid variable range is 0 - 3\n" );
       return 0;
     }
     if (!obj) return( 0 );

     dims[0] = obj->nCell;
     dims[1] = obj->nLev;
     dims[2] = obj->nRow;
     dims[3] = obj->nCol;
     out = (PyArrayObject *) PyArray_FromDims( 4, dims, PyArray_FLOAT64 );
     memcpy( (float64 *) out->data, obj->val0,
	     obj->nCell * obj->cellSize * sizeof( float64 ) );

     return( (PyObject *) out );
  }
  %apply (FMField *in) {
      (FMField *out),
      (FMField *in)
  };
  %apply (int32 *array, int32 len) {
      (int32 *elList, int32 elList_nRow)
  };
  /*!
    @par Revision history:
    - 24.04.2007, c
  */
  int32 integrate( FMField *out, FMField *in, int32 mode=0 ) {
    return( sg_integrate( self, out, in, mode ) );
  }
  /*!
    @par Revision history:
    - 01.11.2007, c
  */
  int32 integrate_chunk( FMField *out, FMField *in,
			 int32 *elList, int32 elList_nRow, int32 mode ) {
    return( sg_integrateChunk( self, out, in, elList, elList_nRow, mode ) );
  }

  /*!
    @par Revision history:
    - 04.05.2007, c
  */
  int32 evaluate_bfbgm( FMField *bfBGR, FMField *ebfBGR,
			float64 *coorIn, int32 nNod, int32 dim,
			int32 *fis, int32 nFa, int32 nFP,
			int32 *conn, int32 nEl, int32 nEP ) {
    return( sg_evaluateBFBGM( self, bfBGR, ebfBGR, coorIn, nNod, dim,
			      fis, nFa, nFP, conn, nEl, nEP ) );
  }

}

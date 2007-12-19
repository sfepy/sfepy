/*!
  @par Revision history:
  - 20.12.2005, c
*/
#ifndef _TERMS_H_
#define _TERMS_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"

void debug_printConn( int32 *conn, int32 nEP );

int32 ele_extractNodalValuesNBN( FMField *out, FMField *in,
				 int32 *conn );

int32 ele_extractNodalValuesDBD( FMField *out, FMField *in,
				 int32 *conn );

END_C_DECLS

#endif /* Header */

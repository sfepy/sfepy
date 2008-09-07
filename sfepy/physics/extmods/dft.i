/* -*- C -*- */
%module dft

%{
#include "f2c.h"
%}

typedef int integer;
typedef float real;

%include "typemaps.i"

%apply (real *OUTPUT) { (real *n) };
%apply (real *INPUT) { (real *vxc) };
%apply (integer *INPUT) { (integer *relat) };

%inline %{
void getvxc(real *n, real *vxc, integer *relat) {
  getvxc_(n, vxc, relat);
}
%}

/* -*- C -*- */
%module dft

typedef int integer;
typedef float real;

%include "typemaps.i"

%apply (real *OUTPUT) { (real *n) };
%apply (real *INPUT) { (real *vxc) };
%apply (integer *INPUT) { (integer *relat) };

%inline %{
void getvxc(double *n, double *vxc, int *relat) {
  getvxc_(*n, vxc, *relat);
}
%}

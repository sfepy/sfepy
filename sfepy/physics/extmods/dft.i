/* -*- C -*- */
%module dft

%{
#include "dft.h"
%}

%rename( "getvxc" ) vxc;
double vxc(double n, int mode, int relat);

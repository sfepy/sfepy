import os.path as op
import re

import numpy as nm

from sfepy.base.ioutils import write_bb

##
# 06.04.2005, c 
def read_vp_fluent( file ):
    desc = file.readline().split()
    ii = desc.index( 'cellnumber' )
    ivx = desc.index( 'x-velocity' )
    ivy = desc.index( 'y-velocity' )
    ivz = desc.index( 'z-velocity' )
    ip = desc.index( 'pressure' )
    print desc
##    print ii, ivx, ivy, ivz, ip

    lines = file.readlines()

    n_row = len( lines )

    vec_v = nm.zeros( (n_row, 3), nm.float64 )
    vec_p = nm.zeros( (n_row, 1), nm.float64 )

    for line in lines:
        fields = re.sub( r'(\d\d)([+-])', r'\1 \2', line ).split()

        ir = int( fields[ii] ) - 1
        vec_v[ir,0] = float( fields[ivx] )
        vec_v[ir,1] = float( fields[ivy] )
        vec_v[ir,2] = float( fields[ivz] )
        vec_p[ir,0] = float( fields[ip] )

    return( (vec_v, vec_p) )

##
# 12.12.2005, c
# 20.01.2006
# 26.07.2006
def read_state( conf, fields, write_meshes = False ):
    # Read state.
    try:
        file = open( conf.filename_vp, "r" );
    except:
        print "Cannot open " + conf.filename_vp + " for reading!";
        raise "ERR_FileOpen"
    vec_v_in, vec_p_in = read_vp_fluent( file )
    file.close()

    print vec_v_in.shape
    print vec_p_in.shape

    field = fields['3_velocity']
    vec_v = field.interp_c_vals_to_n_vals( vec_v_in )
    vec_p = fields['pressure'].interp_c_vals_to_n_vals( vec_p_in )
#    vec_p = None

    if write_meshes:
        trunk = op.splitext( conf.filename_mesh )[0]

        fd = open( trunk + '_3_velocity.bb', 'w' )
        write_bb( fd, vec_v, 2 )
        fd.close()

        fd = open( trunk + '_3_velocity_orig.bb', 'w' )
        write_bb( fd, vec_v_in, 1 )
        fd.close()

        field.write_mesh( trunk + '_%s.mesh' )

    return vec_v, vec_p

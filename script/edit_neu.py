#!/usr/bin/python
# 20.12.2005, c
import sys
"""
Replaces vertices in the Gambit .neu file by the vertices
given in Medit .mesh file.
"""

usage = """
Three args required! (file_in|-, file_mesh, file_out|-)
"""

if (len( sys.argv ) == 4):
    filename_in = sys.argv[1];
    filename_mesh = sys.argv[2];
    filename_out = sys.argv[3];
else:
    print usage
    raise '123'

##
# Read vertices.
fd = open( filename_mesh, "r" ); 
while 1:
    line = fd.readline().split()
    if not len( line ):
        break
    elif (line[0] == 'Vertices'):
        n_nod = int( fd.readline() )
        coors = []
        for ii in range( n_nod ):
            line = [float( ic ) for ic in fd.readline().split()]
            coors.append( '%10d   % 14.10e   % 14.10e   % 14.10e\n'\
                          % ((ii+1,) +  tuple( line[:-1] )) )
fd.close()

if (filename_in == '-'):
    fdi = sys.stdin
else:
    fdi = open( filename_in, "r" ); 

if (filename_out == '-'):
    fdo = sys.stdout
else:
    fdo = open( filename_out, "w" ); 

while 1:
    line = fdi.readline()
    if not len( line ): break
    row = line.split()
    if len( row ) == 0: continue
    if (row[0] == 'NUMNP'):
        fdo.write( line )
        line = fdi.readline()
        row = line.split()
        if n_nod != int( row[0] ):
            raise RuntimeError, 'Meshes not compatible! (%d == %s)'\
                  % (n_nod, row[0])
        fdo.write( line )
    elif (row[0] == 'NODAL'):
        fdo.write( line )
        for ii in range( n_nod ):
            aux = fdi.readline()
            fdo.write( coors[ii] )
    else:
        fdo.write( line )
        
if (filename_in != '-'):
    fdi.close()
if (filename_out != '-'):
    fdo.close()

#!/usr/bin/python
# 20.12.2005, c
import sys
"""
Replaces vertices in the Gambit .neu file by the vertices
given in Medit .mesh file.
"""

usage = """
Three args required! (fileIn|-, fileMesh, fileOut|-)
"""

if (len( sys.argv ) == 4):
    fileNameIn = sys.argv[1];
    fileNameMesh = sys.argv[2];
    fileNameOut = sys.argv[3];
else:
    print usage
    raise '123'

##
# Read vertices.
fd = open( fileNameMesh, "r" ); 
while 1:
    line = fd.readline().split()
    if not len( line ):
        break
    elif (line[0] == 'Vertices'):
        nNod = int( fd.readline() )
        coors = []
        for ii in range( nNod ):
            line = [float( ic ) for ic in fd.readline().split()]
            coors.append( '%10d   % 14.10e   % 14.10e   % 14.10e\n'\
                          % ((ii+1,) +  tuple( line[:-1] )) )
fd.close()

if (fileNameIn == '-'):
    fdi = sys.stdin
else:
    fdi = open( fileNameIn, "r" ); 

if (fileNameOut == '-'):
    fdo = sys.stdout
else:
    fdo = open( fileNameOut, "w" ); 

while 1:
    line = fdi.readline()
    if not len( line ): break
    row = line.split()
    if len( row ) == 0: continue
    if (row[0] == 'NUMNP'):
        fdo.write( line )
        line = fdi.readline()
        row = line.split()
        if nNod != int( row[0] ):
            raise RuntimeError, 'Meshes not compatible! (%d == %s)'\
                  % (nNod, row[0])
        fdo.write( line )
    elif (row[0] == 'NODAL'):
        fdo.write( line )
        for ii in range( nNod ):
            aux = fdi.readline()
            fdo.write( coors[ii] )
    else:
        fdo.write( line )
        
if (fileNameIn != '-'):
    fdi.close()
if (fileNameOut != '-'):
    fdo.close()

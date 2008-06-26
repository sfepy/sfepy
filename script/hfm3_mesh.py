#!/usr/bin/python
# 29.06.2004, c
import sys
import fileinput

if (len( sys.argv ) == 3):
    fileNameIn = sys.argv[1];
    fileNameOut = sys.argv[2];
else:
    print 'Two args required!'
    raise '123'

if (fileNameOut == '-'):
    fileOut = sys.stdout
else:
    fileOut = open( fileNameOut, "w" ); 

fileOut.write( """MeshVersionFormatted 1
# Mesh converted by hfm3_mesh.py

Dimension 3
""" )

mode = 0
nod = []
el = {'brick8' : [], 'tetra4' : []}
for line in fileinput.input( fileNameIn ):
    row = [ii for ii in line.split()]
    if len( row ) == 0: continue
    if (row[0] == 'node3'):
        nod.append( row[1:5] )
    elif (row[0] == 'brick8'):
        el['brick8'].append( row[1:11] )
    elif (row[0] == 'tetra4'):
        el['tetra4'].append( row[1:7] )

fileOut.write( "Vertices\n%d\n" % len( nod ) )
for nn in nod:
    fileOut.write( " ".join( nn[1:4] ) + " 0\n" )

fileOut.write( "Hexahedra\n%d\n" % len( el['brick8'] ) )
for ee in el['brick8']:
    fileOut.write( " ".join( ee[2:11] ) + " " + ee[1] + "\n" )

fileOut.write( "Tetrahedra\n%d\n" % len( el['tetra4'] ) )
for ee in el['tetra4']:
    fileOut.write( " ".join( ee[2:7] ) + " " + ee[1] + "\n" )

if (fileNameOut != '-'):
    fileOut.close()

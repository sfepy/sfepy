#!/usr/bin/python
# 29.06.2004, c
import sys
import fileinput

if (len( sys.argv ) == 3):
    file_name_in = sys.argv[1];
    file_name_out = sys.argv[2];
else:
    print 'Two args required!'
    raise '123'

if (file_name_out == '-'):
    file_out = sys.stdout
else:
    file_out = open( file_name_out, "w" ); 

file_out.write( """MeshVersionFormatted 1
# Mesh converted by hfm3_mesh.py

Dimension 3
""" )

mode = 0
nod = []
el = {'brick8' : [], 'tetra4' : []}
for line in fileinput.input( file_name_in ):
    row = [ii for ii in line.split()]
    if len( row ) == 0: continue
    if (row[0] == 'node3'):
        nod.append( row[1:5] )
    elif (row[0] == 'brick8'):
        el['brick8'].append( row[1:11] )
    elif (row[0] == 'tetra4'):
        el['tetra4'].append( row[1:7] )

file_out.write( "Vertices\n%d\n" % len( nod ) )
for nn in nod:
    file_out.write( " ".join( nn[1:4] ) + " 0\n" )

file_out.write( "Hexahedra\n%d\n" % len( el['brick8'] ) )
for ee in el['brick8']:
    file_out.write( " ".join( ee[2:11] ) + " " + ee[1] + "\n" )

file_out.write( "Tetrahedra\n%d\n" % len( el['tetra4'] ) )
for ee in el['tetra4']:
    file_out.write( " ".join( ee[2:7] ) + " " + ee[1] + "\n" )

if (file_name_out != '-'):
    file_out.close()

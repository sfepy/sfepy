name = '2_3'
vCoors = {'m0' : ( 0, 0 ),
          'm1' : ( 1, 0 ),
          'm2' : ( 0, 1 )}
vEdges = (('m0', 'm1'),
          ('m1', 'm2'),
          ('m2', 'm0'))
vFaces = (('m0', 'm1', 'm2'))
sCoors = {'s0' : ( 0,),
          's1' : ( 1,)}
sEdges = {'s2' : (('s0', 's1'),)}

interpolation = '2_3_P1'

# (root1, vertices of direction vectors, swap from, swap to,
# root2, ...)
orientation = {
     'vVecs' : (0, (1, 2), 1, 2),
}

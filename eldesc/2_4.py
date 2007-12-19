name = '2_4'
vCoors = {'m0' : (-1,-1),
          'm1' : ( 1,-1),
          'm2' : ( 1, 1),
          'm3' : (-1, 1)}
vEdges = (('m0', 'm1'),
          ('m1', 'm2'),
          ('m2', 'm3'),
          ('m3', 'm0'))
vFaces = (('m0', 'm1', 'm2', 'm3'))
sCoors = {'s0' : ( -1,),
          's1' : (  1,)}
sEdges = {'s2' : (('s0', 's1'),)}

interpolation = '2_4_Q1'

# Not finished...
orientation = {
    'vVecs' : (0, (1, 3), (0, 1), (3, 2) ),
}

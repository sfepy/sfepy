name = '3_4'
vCoors = {'m0' : ( 0, 0, 0),
          'm1' : ( 1, 0, 0),
          'm2' : ( 0, 1, 0),
          'm3' : ( 0, 0, 1)}
vEdges = (('m0', 'm1'),
          ('m1', 'm2'),
          ('m2', 'm0'),
          ('m0', 'm3'),
          ('m1', 'm3'),
          ('m2', 'm3'))
vFaces = (('m0', 'm2', 'm1'),
          ('m0', 'm3', 'm2'),
          ('m0', 'm1', 'm3'),
          ('m1', 'm2', 'm3'))
sCoors = {'s0' : ( 0, 0),
          's1' : ( 1, 0),
          's2' : ( 0, 1)}
sEdges = {'s3' : (('s0', 's1'),
                  ('s1', 's2'),
                  ('s2', 's0'))}
sFaces = {'s3' : ('s0', 's1', 's2')}

interpolation = '3_4_P1'

# (root1, vertices of direction vectors, swap from, swap to,
# root2, ...)
orientation = {
     'vVecs' : (0, (1, 2, 3), 0, 3),
}

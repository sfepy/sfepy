name = '3_4'
v_coors = {'m0' : ( 0, 0, 0),
          'm1' : ( 1, 0, 0),
          'm2' : ( 0, 1, 0),
          'm3' : ( 0, 0, 1)}
v_edges = (('m0', 'm1'),
          ('m1', 'm2'),
          ('m2', 'm0'),
          ('m0', 'm3'),
          ('m1', 'm3'),
          ('m2', 'm3'))
v_faces = (('m0', 'm2', 'm1'),
          ('m0', 'm3', 'm2'),
          ('m0', 'm1', 'm3'),
          ('m1', 'm2', 'm3'))
s_coors = {'s0' : ( 0, 0),
          's1' : ( 1, 0),
          's2' : ( 0, 1)}
s_edges = {'s3' : (('s0', 's1'),
                  ('s1', 's2'),
                  ('s2', 's0'))}
s_faces = {'s3' : ('s0', 's1', 's2')}

interpolation = '3_4_P1'

# (root1, vertices of direction vectors, swap from, swap to,
# root2, ...)
orientation = {
     'v_vecs' : (0, (1, 2, 3), 0, 3),
}

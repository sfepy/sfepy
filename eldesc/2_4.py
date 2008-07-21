name = '2_4'
v_coors = {'m0' : (-1,-1),
          'm1' : ( 1,-1),
          'm2' : ( 1, 1),
          'm3' : (-1, 1)}
v_edges = (('m0', 'm1'),
          ('m1', 'm2'),
          ('m2', 'm3'),
          ('m3', 'm0'))
v_faces = (('m0', 'm1', 'm2', 'm3'))
s_coors = {'s0' : ( -1,),
          's1' : (  1,)}
s_edges = {'s2' : (('s0', 's1'),)}

interpolation = '2_4_Q1'

# Not finished...
orientation = {
    'v_vecs' : (0, (1, 3), (0, 1), (3, 2) ),
}

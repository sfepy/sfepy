name = '2_3'
v_coors = {'m0' : ( 0, 0 ),
          'm1' : ( 1, 0 ),
          'm2' : ( 0, 1 )}
v_edges = (('m0', 'm1'),
          ('m1', 'm2'),
          ('m2', 'm0'))
v_faces = (('m0', 'm1', 'm2'))
s_coors = {'s0' : ( 0,),
          's1' : ( 1,)}
s_edges = {'s2' : (('s0', 's1'),)}

interpolation = '2_3_P1'

# (root1, vertices of direction vectors, swap from, swap to,
# root2, ...)
orientation = {
     'v_vecs' : (0, (1, 2), 1, 2),
}

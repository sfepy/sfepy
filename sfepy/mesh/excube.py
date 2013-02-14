import math
import numpy as np

class mesh_writer:
    def __init__(self, filename, dim = 3):
        if not isinstance(filename, file):
           self.file = open(filename,"w")
           self.close = True
        else:
           self.file = filename
           self.close = False
        self.dim = 3
        self.elementType = 'Hexahedra'
        self.verticles = 8
        self.nodes_cnt = [0]
        self.elements_cnt = [0]
        self.group = 0
        self.nodes = [np.empty((1000, self.dim),np.double)]
        self.elements = [np.empty((1000, self.verticles),np.integer)]
#        self.debug=True
        self.output=True
        self.debug=False

    def node(self, *args):
        if(self.nodes_cnt[self.group] + 1 == self.nodes[self.group].shape[0]):
           self.nodes[self.group]=np.tile(self.nodes[self.group], (3,1))
        self.nodes[self.group][self.nodes_cnt[self.group],:]=args
        self.nodes_cnt[self.group]+=1

    def printnode(self, i, group=None):
        if group is None: group=self.group
        print self.nodes[group][i]

    def printelement(self, i, group=None):
        if group is None: group=self.group
        print self.nodes[group][self.elements[group][i]]
        print self.elements[group][i]
 
    def ex_element(self, args, reverse_half=False, reverse_2half=False): 
        if(reverse_half):
           half=len(args)/2
           a=list(args[0:half])
           a.reverse()
           args=a + list(args[half:])
        if(reverse_2half):
           half=len(args)/2
           a=list(args[half:])
           a.reverse()
           args=list(args[:half])+a
        self.element(*args)

    def element(self, *args):
        if(self.elements_cnt[self.group] +1  == self.elements[self.group].shape[0]):
           self.elements=np.tile(self.elements,(3,1))
        self.elements[self.group][self.elements_cnt[self.group],:]=args
        if self.debug:
           self.check(self.elements[self.group][self.elements_cnt[self.group],:])
        self.elements_cnt[self.group]+=1

    def stage(self, stage):
        if self.output:
           g=self.group
           print "STAGE: %s  -  group: %i, nodes: %i, elements: %i" % \
                 (stage, g, self.nodes_cnt[g], self.elements_cnt[g])

    def line(self, a=''):
        self.file.write(str(a) + "\n")

    def write_numbers(self, a):
        for x in a:
          self.file.write(str(x) + ' ')

    def to_stdout(self):
        import sys
        tmp=self.file
        self.file=sys.stdout
        self.write()
        self.file=tmp      

    def write(self):
        self.line("MeshVersionFormated 1")
        self.line("Dimension %i" % self.dim)
        self.line("Vertices")
        self.line(sum(self.nodes_cnt))
        u=1
        for i,nodes in enumerate(self.nodes):
          for x in nodes[0:self.nodes_cnt[i]]:
            self.write_numbers(x)
            self.line(i)

        self.line(self.elementType)
        self.line(sum(self.elements_cnt))
        for i,elements in enumerate(self.elements):
          for x in elements[0:self.elements_cnt[i]]:
              self.write_numbers(x+1)
              self.line(i)
        if self.close: self.file.close()
 


    def check(self, e):
                from numpy.linalg import inv
                import exceptions
                 
                def inTetra(a,b,c,d,p):
                    try:
                      T=inv(np.vstack((a-d,b-d,c-d)).T)
                      l=T.dot((p-d))
                    except Exception, e:
                      print "SINGULAR MATRIX"
                      return True

                    out=(l>=0.0).all() and (l<=1.0).all() and l.sum<=1
                    return out

                def inTet(a,b,c,d,p):
                    nodes=self.nodes[self.group]
                    return inTetra(nodes[a], nodes[b], nodes[c],nodes[d], p)

                def ccheck(e, p):
                    """p=np.array((0.5,0.5,0.5))
                    self.nodes[e[0]]=(0,0,0)
                    self.nodes[e[1]]=(1,0,0)
                    self.nodes[e[2]]=(1,1,0)
                    self.nodes[e[3]]=(0,1,0)
                    self.nodes[e[4]]=(0,0,1)
                    self.nodes[e[5]]=(1,0,1)
                    self.nodes[e[6]]=(1,1,1)
                    self.nodes[e[7]]=(0,1,1) """
                    if (e >= self.nodes_cnt[self.group]).any() or (e<0).any():
                        print "INVALID NODE ", e
                        exit()
                    for x in xrange(0,4):
                        for y in xrange(0,4):
                            if inTet(e[x],e[(x+1) % 4], e[(x+2) % 4], e[y+4], p) or \
                               inTet(e[4+x],e[(x+1) % 4+4], e[(x+2) % 4+4], e[y], p):
                                  print e[0:4]
                                  print e[4:8]
                                  print nodes[e[0:4]]
                                  print nodes[e[4:8]]
                                  exit()
                                  return True
                    #exit()
                    return False
                nodes=self.nodes[self.group]
                mn=nodes[e[0]].copy()
                mx=mn.copy()
                
                for z in e[1:]:
                    z-=1
                    for i in xrange(0,self.dim):
                        if mn[i]>nodes[z,i]: mn[i]=nodes[z,i]
                        if mx[i]<nodes[z,i]: mx[i]=nodes[z,i]
                u=0
                for n in nodes[0:self.nodes_cnt[self.group]]:
                    u+=1
                    if (u==e).any(): continue
                    if (mn<=n).all()  and(mx>=n).all():
                       if ccheck(e-1,n):

                          self.write()
                          print "IN ELEMENT: ",
                          for x in e:
                              print (str(x) + ' '),
                          print "\n WITH BOUNDING BOX %d, %d, %d : %d, %d, %d " % (mn[0],mn[1], mn[2], mx[0],mx[1], mx[2])
                          for x in e:
                              print self.nodes[x-1]
                          print " IS NODE %i \n" % u ,
                          print n
                          exit()


def excube(size, a, b, mesh=None):
    def normalize(what):
        if not isinstance(what, np.ndarray):
           what=np.array(what, np.float)
        if what.size == 1:
           what=np.tile(np.array(what, np.float) , 3)
        return what
    if mesh is None:
       import sys
       mesh=sys.stdout

    if isinstance(mesh, (str, file)):
       mesh=mesh_writer(mesh)

    size = normalize(size)
    a = normalize(a)
    b = normalize(b)

    cnt = np.array(np.floor(a / size), int)
    step = a / cnt


    cube=np.empty((8), np.integer)
    mesh.stage('Inner nodes')
    for z in xrange(-cnt[2], cnt[2]+1):
        for y in xrange(-cnt[1], cnt[1]+1):
            for x in xrange(-cnt[0], cnt[0]+1):
                mesh.node(x*step[0], y*step[1], z*step[2])

    start = (cnt*2 + 1).prod()
    cube[0]=0
    cube[1]=start

    #hyper=(b-a).sum()/3
    #hyper=integer(hyper/size/2)
    #hyper=np.linspace(0,1.0/3,hyper) * np.linspace(0.0,3.0,hyper)

    #prumerna delka
    norm = a.sum() / 3.0
    s_n = step.sum()/ 3.0 / norm
    b_n = b.sum() / 3.0 / norm
    r = b / a

    #solving recurent equation for nodes vertical distance is equal to
    #(average) horizontal distance give this equation for item heigts
    #q = (2 + ratio*s_n - s_n) / (2 - ratio*s_n + s_n)
    #f = lambda n: (1/(1-ratio))  * (q**n) - 1/(1-ratio)
    q = (2 + s_n)/(2 - s_n)
    f = lambda n: q**n
    last=max(1,int(math.log(b_n + 1,q)))
    #print last,math.log(1+(1-ratio),q)
    #last=10
    corr=f(last)


    coefs=np.fromiter((
           (f(x)-1.0)/(corr-1.0)
          for x in xrange(1,last+1)), np.double, last)

    """ EDGE CUBES
            +Z
            1  +Y
            | 3
            |/
        6---0---5+X
           /|
          4 |
            2
       1 & 2 direction 0,0,1
       3 & 4 direction 0,1,0
       5 & 6 direction 1,0,0

       1,3,5 are in positive direction
       2,4,6 negative


    """
    #1,2
    mesh.stage('Outer nodes 1,2')
    for s in [1, -1]:
        for z in coefs:
            for y in xrange(-cnt[1], cnt[1]+1):
                for x in xrange(-cnt[0], cnt[0]+1):
                    mesh.node( x*step[0] * (1+r[0]*z), y*step[1] * (1+r[0]*z), s*(b[2]*z + a[2]) )
#                   print x*step[0] * (1+z), y*step[1] * (1+z), s*(b[2]*z + a[2])
    #print s, b[2], z, a[2]
    nodes=(2*cnt[0]+1) * (2*cnt[1] + 1) * coefs.size
    cube[2]=cube[1] + nodes
    cube[3]=cube[2] + nodes

    #3,4
    mesh.stage('Outer nodes 3,4')
    for s in [1, -1]:
        for z in coefs:
            for y in xrange(-cnt[2]+1, cnt[2]):
                for x in xrange(-cnt[0], cnt[0]+1):
                    mesh.node( x*step[0] * (1+r[0]*z), s*(b[1]*z + a[1]), y*step[2] * (1+r[2]*z)  )
    nodes=(2*cnt[0] + 1) * (2*cnt[2] - 1) * coefs.size
    cube[4]=cube[3] + nodes
    cube[5]=cube[4] + nodes

    #5,6
    mesh.stage('Outer nodes 5,6')
    for s in [1, -1]:
        for z in coefs:
            for y in xrange(-cnt[2]+1, cnt[2]):
                for x in xrange(-cnt[1]+1, cnt[1]):
                    mesh.node( s*(b[0]*z + a[0]), x*step[1] * (1+r[1]*z), y*step[2] * (1+r[2]*z)  )
    nodes=(2*cnt[1] - 1) * (2*cnt[2] - 1) * coefs.size
    cube[6]=cube[5] + nodes
    cube[7]=cube[6] + nodes

    #
    #inner elements in cube
    #
    ff=range(0,7)
    mesh.stage('cube 0')
    #basic cube 0
    bdim = cnt[0]*2+1
    cdim = bdim * (cnt[1]*2+1)
    f = lambda x,y,z, bdim=bdim, cdim=cdim: cnt[0] + x + (y +cnt[1]) * bdim + (z+cnt[2]) * cdim
    ff[0]=f
    for z in xrange(-cnt[2], cnt[2]):
        for y in xrange(-cnt[1], cnt[1]):
            for x in xrange(-cnt[0], cnt[0]):
                mesh.element(
                    f(x,y,z),
                    f(x+1,y,z),
                    f(x+1,y+1,z),
                    f(x,y+1,z),
                    f(x,y,z+1),
                    f(x+1,y,z+1),
                    f(x+1,y+1,z+1),
                    f(x,y+1,z+1)
                  )

    #cube 1,2
    bdim = cnt[0]*2+1
    cdim = bdim * (cnt[1]*2+1)

    for s in (1,2):
        mesh.stage('cube %i' %s)
        f = lambda x,y,z, bdim=bdim, cdim=cdim, s=s:  cnt[0] + x + (y+cnt[1]) * bdim + z * cdim + cube[s]
        ff[s]=f
        for z in xrange(0, coefs.size-1):
            for x in xrange(-cnt[0], cnt[0]):
                for y in xrange(-cnt[1], cnt[1]):
                    mesh.element(
                       f(x,y,z),
                       f(x+1,y,z),
                       f(x+1,y+1,z),
                       f(x,y+1,z),
                       f(x,y,z+1),
                       f(x+1,y,z+1),
                       f(x+1,y+1,z+1),
                       f(x,y+1,z+1)
                     )
    
    #cube 3,4
    bdim = cnt[0]*2+1
    cdim = bdim * (cnt[2]*2-1)

    for s in (3,4):
        mesh.stage('cube %i' %s)
        f = lambda x,y,z, bdim=bdim, cdim=cdim, s=s: cnt[0] + x + (y + cnt[2]-1) * bdim + z * cdim + cube[s]
        ff[s]=f
       
        for z in xrange(0, coefs.size-1):
            for y in xrange(-cnt[2]+1, cnt[2]-1):
                for x in xrange(-cnt[0], cnt[0]):
                    mesh.element(
                       f(x,y,z),
                       f(x+1,y,z),
                       f(x+1,y+1,z),
                       f(x,y+1,z),
                       f(x,y,z+1),
                       f(x+1,y,z+1),
                       f(x+1,y+1,z+1),
                       f(x,y+1,z+1)
                     )

    #cube 5,6
    bdim = cnt[1]*2-1
    cdim = bdim * (cnt[2]*2-1)
    for s in (5,6):
        mesh.stage('cube %i' %s)
        f = lambda x,y,z, bdim=bdim, cdim=cdim, s=s:  x + cnt[2] -1 + (y +cnt[2]-1) * bdim + z * cdim + cube[s]
        ff[s]=f
        for z in xrange(0, coefs.size-1):
            for y in xrange(-cnt[2]+1, cnt[2]-1):
                for x in xrange(-cnt[1]+1, cnt[1]-1):
                    mesh.element(
                       f(x,y,z),
                       f(x+1,y,z),
                       f(x+1,y+1,z),
                       f(x,y+1,z),
                       f(x,y,z+1),
                       f(x+1,y,z+1),
                       f(x+1,y+1,z+1),
                       f(x,y+1,z+1)
                     )
    f=ff


    #face 0:1,2
    for s,z in zip((1,2), (cnt[2], -cnt[2])):
        mesh.stage('face 0:%i' %s)
        for x in xrange(-cnt[0], cnt[0]):
            for y in xrange(-cnt[1], cnt[1]):
                mesh.element(
                       f[0](x,y,z),
                       f[0](x+1,y,z),
                       f[0](x+1,y+1,z),
                       f[0](x,y+1,z),
                       f[s](x,y,0),
                       f[s](x+1,y,0),
                       f[s](x+1,y+1,0),
                       f[s](x,y+1,0)
                       )


   #face 0:3,4
    for s,y in zip((3,4), (cnt[1], -cnt[1])):
        mesh.stage('face 0:%i' %s)
        for x in xrange(-cnt[0], cnt[0]):
            for z in xrange(-cnt[2]+1, cnt[2]-1):
                mesh.element(
                       f[0](x,y,z),
                       f[0](x+1,y,z),
                       f[0](x+1,y,z+1),
                       f[0](x,y,z+1),
                       f[s](x,z,0),
                       f[s](x+1,z,0),
                       f[s](x+1,z+1,0),
                       f[s](x,z+1,0)
                       )

    #face 0:5,6
    for s,x in zip((5,6), (cnt[0], -cnt[0])):
        mesh.stage('face 0:%i' %s)
        for y in xrange(-cnt[1]+1, cnt[1]-1):
            for z in xrange(-cnt[2]+1, cnt[2]-1):
                mesh.element(
                       f[0](x,y,z),
                       f[0](x,y+1,z),
                       f[0](x,y+1,z+1),
                       f[0](x,y,z+1),
                       f[s](y,z,0),
                       f[s](y+1,z,0),
                       f[s](y+1,z+1,0),
                       f[s](y,z+1,0)
                       )



    #
    #inner edges
    #
    #edge (3,4),(1,2),0
    for s,y in zip((3,4), (cnt[1], -cnt[1])):
        for ss,z in zip((1,2), (1, -1)):
            mesh.stage('edge %i:%i' % (s,ss))
            for x in xrange(-cnt[0], cnt[0]):
                mesh.element(
                       f[0](x,y,z*cnt[2]),
                       f[0](x+1,y,z*cnt[2]),
                       f[0](x+1,y,z*(cnt[2]-1)),
                       f[0](x,y,z*(cnt[2]-1)),

                       f[ss](x,y, 0),
                       f[ss](x+1,y, 0),
                       f[s](x+1,z*(cnt[2]-1), 0),
                       f[s](x,z*(cnt[2]-1), 0),
                     )
    #
    #inner edges
    #
    #edge (5,6),(1,2),0
    for s,x in zip((5,6), (cnt[0], -cnt[0])):
        for ss,z in zip((1,2), (1, -1)):
            mesh.stage('edge %i:%i' % (s,ss))
            for y in xrange(-cnt[1]+1, cnt[1]-1):
                mesh.element(
                       f[0](x,y,z*cnt[2]),
                       f[0](x,y+1,z*cnt[2]),
                       f[0](x,y+1,z*(cnt[2]-1)),
                       f[0](x,y,z*(cnt[2]-1)),

                       f[ss](x,y, 0),
                       f[ss](x,y+1, 0),
                       f[s](y+1,z*(cnt[2]-1), 0),
                       f[s](y,z*(cnt[2]-1), 0),
                     )

    #
    #inner edges
    #
    #edge (5,6),(3,4),0
    for s,x in zip((5,6), (cnt[0], -cnt[0])):
        for ss,y in zip((3,4), (1, -1)):
            mesh.stage('edge %i:%i' % (s,ss))
            for z in xrange(-cnt[2]+1, cnt[2]-1):
                mesh.element(
                       f[0](x,y*cnt[1],z),
                       f[0](x,y*cnt[1],z+1),
                       f[0](x,y*(cnt[1]-1),z+1),
                       f[0](x,y*(cnt[1]-1),z),

                       f[ss](x,z, 0),
                       f[ss](x,z+1, 0),
                       f[s](y*(cnt[1]-1),z+1, 0),
                       f[s](y*(cnt[1]-1),z, 0),
                     )



    #exit()
    #
    #verticles
    #
    #0,0,0,0,(1,2),(1,2),(3,4),(5,6)
    for ss,z in zip((1,2), (1, -1)):
       for s,y in zip((3,4), (1, -1)):
          for sss,x in zip((5,6), (cnt[0], -cnt[0])):
                mesh.stage('verticle %i:%i:%i' % (s,ss,sss))
                mesh.element(
                       f[0](x,y*cnt[1],z*cnt[2]),
                       f[0](x,y*(cnt[1]-1),z*cnt[2]),
                       f[0](x,y*(cnt[1]-1),z*(cnt[2]-1)),
                       f[0](x,y*cnt[1],z*(cnt[2]-1)),

                       f[ss](x,y*cnt[1], 0),
                       f[ss](x,y*(cnt[1]-1), 0),
                       f[sss](y*(cnt[1]-1),z*(cnt[2]-1), 0),
                       f[s](x,z*(cnt[2]-1), 0),
                     )


    #
    #outer faces
    #
    #
    #(3,4),(1,2)
    for s, bz in zip((1,2), (cnt[2]-1, -cnt[2]+1)):
       for ss, ay in zip((3,4), (cnt[1], -cnt[1])):
          mesh.stage('face %i:%i' % (s,ss))
          for z in xrange(0,coefs.size-1):
              for x in xrange(-cnt[0], cnt[0]):
                mesh.element(
                       f[s](x,ay,z),
                       f[s](x,ay,z+1),
                       f[s](x+1,ay,z+1),
                       f[s](x+1,ay,z),

                       f[ss](x,bz,z),
                       f[ss](x,bz,z+1),
                       f[ss](x+1,bz,z+1),
                       f[ss](x+1,bz,z),
                     )

    #(5,6),(1,2)
    for s, bz in zip((1,2), (cnt[2]-1, -cnt[2]+1)):
       for ss, ax in zip((5,6), (cnt[0], -cnt[0]) ):
          mesh.stage('face %i:%i' % (s,ss))
          for y in xrange(-cnt[1]+1, cnt[1]-1):
             for z in xrange(0,coefs.size-1):
                mesh.element(
                       f[s](ax,y,z),
                       f[s](ax,y+1,z),
                       f[s](ax,y+1,z+1),
                       f[s](ax,y,z+1),

                       f[ss](y,bz,z),
                       f[ss](y+1,bz,z),
                       f[ss](y+1,bz,z+1),
                       f[ss](y,bz,z+1),
                     )

    #(5,6),(3,4)
    for s, by in zip((3,4), (cnt[1]-1, -cnt[1]+1)):
       for ss, ax in zip((5,6), (cnt[0], -cnt[0]) ):
          mesh.stage('face %i:%i' % (s,ss))
          for x in xrange(-cnt[2]+1, cnt[2]-1):
             for z in xrange(0,coefs.size-1):
                mesh.element(
                       f[s](ax,x,z),
                       f[s](ax,x,z+1),
                       f[s](ax,x+1,z+1),
                       f[s](ax,x+1,z),

                       f[ss](by,x,z),
                       f[ss](by,x,z+1),
                       f[ss](by,x+1,z+1),
                       f[ss](by,x+1,z),
                     )

    #
    # Outer edges
    #
    #1,1,1,1,3,3,5,5
    for s, az in zip((1,2), (cnt[2]-1, -cnt[2]+1)):
       for ss, by in zip((3,4), (1,-1) ):
          for sss, cz in zip((5,6), (1, -1)):
             mesh.stage('edge %i:%i:%i' % (s,ss,sss))
             for z in xrange(0,coefs.size-1):
                mesh.element(
                       f[s](cz*cnt[0],by*cnt[1] ,z),
                       f[s](cz*cnt[0],by*cnt[1],z+1),
                       f[s](cz*cnt[0],by*(cnt[1]-1),z+1),
                       f[s](cz*cnt[0],by*(cnt[1]-1),z),

                       f[ss](cz*cnt[0],az,z),
                       f[ss](cz*cnt[0],az,z+1),
                       f[sss](by*(cnt[1]-1),az,z+1),
                       f[sss](by*(cnt[1]-1),az,z),
                     )
#    mesh.printelement(5002)
    mesh.write()

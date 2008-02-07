class geometry(object):
    """The geometry is given by a sets of points (d0), lines (d1), surfaces
    (d2) and volumes (d3). A lines are constructed from 2 points, a surface from
    any number of lines, a volume from any number of surfaces.

    Physical volumes are contruted from any number of volumes.


    The self.d0, self.d1, self.d2 and self.d3 are dictionaries holding a map

    geometry element number  ->   instance of point,line,surface of volume

    Example:
    ========

    To get all the points which define a surface 5, use:

    self.d2[5].getpoints()

    This would give you a list [..] of point() instances.
    """
    def __init__(self):
        self.d0={}
        self.d1={}
        self.d2={}
        self.d3={}
        self.phys2={}
        self.phys3={}
    def addpoint(self,n,p):
        "p=[x,y,z]"
        o=point(self,n,p)
        self.d0[o.getn()]=o
    def addline(self,n,l):
        "l=[p1,p2]"
        o=line(self,n,l)
        self.d1[o.getn()]=o
    def addsurface(self,n,s):
        "s=[l1,l2,l3,...]"
        o=surface(self,n,s)
        self.d2[o.getn()]=o
    def addvolume(self,n,v):
        "v=[s1,s2,s3,...]"
        o=volume(self,n,v)
        self.d3[o.getn()]=o
    def addphysicalsurface(self,n,surfacelist):
        "surfacelist=[v1,v2,v3,...]"
        o=physicalsurface(self,n,surfacelist)
        self.phys2[o.getn()]=o
    def addphysicalvolume(self,n,volumelist):
        "volumelist=[v1,v2,v3,...]"
        o=physicalvolume(self,n,volumelist)
        self.phys3[o.getn()]=o
    def getBCnum(self,snum):
        for x in self.phys2:
            if snum in self.phys2[x].surfaces:
                return x
        return 0
    def printinfo(self):
        print "General geometry information:"
        print "  points:",len(self.d0)
        print "  lines:",len(self.d1)
        print "  surfaces:",len(self.d2)
        print "  volumes:",len(self.d3)
        print "Physical entities:"
        print "  surfaces (boundary conditions):"
        for d in self.phys2.values():
            print "    %d: surface numbers %r"%(d.getn(),d.surfaces)
        print "  volumes (regions):"
        for d in self.phys3.values():
            print "    %d: volume numbers %r"%(d.getn(),d.volumes)
    def leaveonlyphysicalvolumes(self):
        points={}
        lines={}
        surfaces={}
        volumes={}
        for e in self.phys3:
            for v in self.phys3[e].getvolumes():
                volumes[v.getn()]=v
                for s in v.getsurfaces():
                    surfaces[s.getn()]=s
                    for l in s.getlines():
                        lines[l.getn()]=l
                        for p in l.getpoints():
                            points[p.getn()]=p
        self.d0=points
        self.d1=lines
        self.d2=surfaces
        self.d3=volumes

class geomobject(object):
    def getn(self):
        return self.n

class point(geomobject):
    def __init__(self,g,n,p):
        self.geom=g
        self.n=n
        self.p=p
    def __add__(self,p):
        return point(self.geom,-1,[a+b for a,b in zip(self.p,p.p)])
    def __sub__(self,p):
        return point(self.geom,-1,[a-b for a,b in zip(self.p,p.p)])
    def __div__(self,num):
        return point(self.geom,-1,[a/num for a in self.p])
    def __mul__(self,num):
        return point(self.geom,-1,[a*num for a in self.p])
    def __rmul__(self,num):
        return self.__mul__(num)
    def getxyz(self):
        return self.p
    def getstr(self):
        return "%f, %f, %f"%self.getxyz()

class line(geomobject):
    def __init__(self,g,n,l):
        self.geom=g
        self.n=n
        self.points=l
    def getpoints(self):
        return [self.geom.d0[x] for x in self.points]

class surface(geomobject):
    def __init__(self,g,n,s):
        self.geom=g
        self.n=n
        self.lines,self.holes=self.separate(s)
    def separate(self,s):
        #FIXME - this is just a quick hack to satisfy all the examples
        if len(s)<=4:
            return s,[]
        elif len(s)==8:
            return s[:4],[s[4:]]
        else:
            return s,[]
    def getlines(self):
        return [self.geom.d1[abs(x)] for x in self.lines]
    def getpoints(self):
        #self.lines contains the numbers of all the lines
        def p(idx):
            "Return the correct point of the line 'idx'"
            if idx>0:
                return self.geom.d1[idx].getpoints()[0]
            else:
                return self.geom.d1[-idx].getpoints()[1]
        return [p(x) for x in self.lines]
    def getholepoints(self):
        def p(idx):
            "Return the correct point of the line 'idx'"
            if idx>0:
                return self.geom.d1[idx].getpoints()[0]
            else:
                return self.geom.d1[-idx].getpoints()[1]
        r=[]
        for hole in self.holes:
            r.append([p(x) for x in hole])
        return r
    def getinsidepoint(self):
        pts=self.getpoints()[:3]
        return (pts[0]+pts[1]+pts[2])/3

class volume(geomobject):
    def __init__(self,g,n,v):
        self.geom=g
        self.n=n
        self.surfaces=v
    def getsurfaces(self):
        return [self.geom.d2[abs(x)] for x in self.surfaces]
    def getinsidepoint(self):
        sfs=self.getsurfaces()[:3]
        pts=[s.getinsidepoint() for s in sfs]
        p0=sfs[0].getpoints()[0]
        direct=(pts[0]+pts[1]+pts[2])/3-p0
        return p0+0.001*direct

class physicalsurface(geomobject):
    def __init__(self,g,n,s):
        self.geom=g
        self.n=n
        self.surfaces=s
    def getsurfaces(self):
        return [self.geom.d2[x] for x in self.surfaces]

class physicalvolume(geomobject):
    def __init__(self,g,n,v):
        self.geom=g
        self.n=n
        self.volumes=v
    def getvolumes(self):
        return [self.geom.d3[x] for x in self.volumes]

from __future__ import print_function
from __future__ import absolute_import
import numpy as nm
import six
from six.moves import range

class geometry(object):
    """The geometry is given by a sets of points (d0), lines (d1), surfaces
    (d2) and volumes (d3). A lines are constructed from 2 points, a surface from
    any number of lines, a volume from any number of surfaces.

    Physical volumes are contruted from any number of volumes.


    The self.d0, self.d1, self.d2 and self.d3 are dictionaries holding a map

    geometry element number  ->   instance of point,line,surface of volume

    Examples
    --------

    To get all the points which define a surface 5, use:

    self.d2[5].getpoints()

    This would give you a list [..] of point() instances.
    """
    def __init__(self, dim=3):
        self.dim = dim
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
    def addpoints(self,ps,off=1):
        "ps=[p1, p2, ...]"
        for i, p in enumerate(ps):
            self.addpoint(i + off, p)
    def addline(self,n,l):
        "l=[p1,p2]"
        o=line(self,n,l)
        self.d1[o.getn()]=o
    def addlines(self,ls,off=1):
        "ls=[l1, l2, ...]"
        for i, l in enumerate(ls):
            self.addline(i + off, l)
    def addsurface(self,n,s, is_hole=False):
        "s=[l1,l2,l3,...]"
        o=surface(self,n,s, is_hole)
        self.d2[o.getn()]=o
    def addsurfaces(self,ss,off=1):
        "s=[s1,s2,s3,...]"
        for i, s in enumerate(ss):
            self.addsurface(i + off, s)
    def addvolume(self,n,v):
        "v=[s1,s2,s3,...]"
        o=volume(self,n,v)
        self.d3[o.getn()]=o
    def addvolumes(self,vs,off=1):
        "v=[v1,v2,v3,...]"
        for i, v in enumerate(vs):
            self.addvolume(i + off, v)
    def addphysicalsurface(self,n,surfacelist):
        "surfacelist=[s1,s2,s3,...]"
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
    def printinfo(self, verbose=False):
        print("General geometry information:")
        print("  dimension:", self.dim)
        print("  points:", len(self.d0))
        if verbose:
            for k, v in six.iteritems(self.d0):
                print("    %d - %s" % (k, v.getstr()))
        print("  lines:", len(self.d1))
        if verbose:
            for k, v in six.iteritems(self.d1):
                print("    %d - " % k, v.points)
        print("  surfaces:", len(self.d2))
        if verbose:
            for k, v in six.iteritems(self.d2):
                if v.is_hole:
                    aux = '(hole)'
                else:
                    aux = ''
                print("    %d%s - " % (k, aux), v.lines)
        print("  volumes:", len(self.d3))
        if verbose:
            for k, v in six.iteritems(self.d3):
                print("    %d - " % k, v.surfaces)
        print("Physical entities:")
        if self.dim == 2:
            print("  surfaces (regions):")
            for d in self.phys2.values():
                print("    %d: surface numbers %r"%(d.getn(),d.surfaces))
        elif self.dim == 3:
            print("  surfaces (boundary conditions):")
            for d in self.phys2.values():
                print("    %d: surface numbers %r"%(d.getn(),d.surfaces))
            print("  volumes (regions):")
            for d in self.phys3.values():
                print("    %d: volume numbers %r"%(d.getn(),d.volumes))

    def leaveonlyphysicalsurfaces(self):
        points={}
        lines={}
        surfaces={}
        volumes={}
        for e in self.phys2:
            for s in self.phys2[e].getsurfaces():
                surfaces[s.getn()]=s
                for l in s.getlines():
                    lines[l.getn()]=l
                    for p in l.getpoints():
                        points[p.getn()]=p
        self.d0=points
        self.d1=lines
        self.d2=surfaces
        self.d3=volumes

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

    def splitlines(self, ls, n):
        repl = {}
        for il in ls:
            l = self.d1[il]
            pts = l.getpoints()
            dp = pts[1] - pts[0]
            t = nm.linspace(0, 1, n + 1)
            points = [pts[0].n, pts[1].n]
            for ii, it in enumerate(t[1:-1]):
                pid = il * 1000 + ii
                self.addpoint(pid, (pts[0] + dp * it).getxyz())
                points.insert(-1, pid)

            lines = []
            for ii in range(n):
                lid = il * 1000 + ii
                self.addline(lid, [points[ii], points[ii + 1]])
                lines.append(lid)

            for s in six.itervalues(self.d2):
                try:
                    idx = s.lines.index(l.n)
                except ValueError:
                    continue

                s.lines.pop(idx)
                for ii, j in enumerate(lines):
                    s.lines.insert(idx + ii, j)

            repl[l.n] = lines
            self.d1.pop(l.n)

    def to_poly_file(self, filename):
        """
        Export geometry to poly format (tetgen and triangle geometry format).

        Parameters
        ----------
        geo : geometry
            geometry description
        filename : string
            file name
        """

        def getinsidepoint(pts):
            direct = (pts[0] + pts[1] + pts[2]) / 3 - pts[0]
            return pts[0] + 0.001 * direct

        if self.dim == 2:
            self.leaveonlyphysicalsurfaces()
        if self.dim == 3:
            self.leaveonlyphysicalvolumes()

        # write nodes
        nodes = []
        map = {}
        for x in self.d0.values():
            assert isinstance(x, point)
            nodes.append(x.getxyz())
            map[x.getn()] = len(nodes)


        s = "# nodes\n%d %d 0 0\n" % (len(nodes), self.dim)
        if self.dim == 2:
            ptstr = " %d %f %f\n"
            ptstr2 = " %d %f %f %d\n"
        else:
            ptstr = " %d %f %f %f\n"
            ptstr2 = " %d %f %f %f %d\n"

        for n, x in enumerate(nodes):
            s += ptstr % tuple([n + 1] + list(x[:self.dim]))

        # facets
        # first write external polygon, then hole polygons and then point in each
        # hole polygon
        facets = []
        if self.dim == 2:

            hole_pts = []
            regions=[]
            for x2 in self.d2.values():
                assert isinstance(x2, surface)
                for x1 in x2.getlines():
                    assert isinstance(x1, line)
                    p = [map[y.getn()] for y in x1.getpoints()]
                    bc = self.getBCnum(x1.getn())
                    facets.append((p, bc))

                for hole in x2.getholepoints():
                    hole_pts.append(hole.getxyz())

        # regions
        for x in self.phys2.values():
            assert isinstance(x, physicalsurface)
            for x2 in x.getsurfaces():
                if not x2.is_hole:
                    regions.append(x2.getinsidepoint().getxyz() + [x.getn()])

            # number of facets, boundary markers=yes
            s += "# segments\n%d 1\n" % len(facets)
            for ii, (p, bc) in enumerate(facets):
                # number of corners, corner 1, corner 2, ...
                s += " %d %s %d\n" % (ii + 1, ' '.join([str(ii) for ii in p]), bc)
            # holes
            s += "# holes\n%d\n" % len(hole_pts)
            for ii, x0 in enumerate(hole_pts):
                # number of corners, corner 1, corner 2, ...
                s += " %d %s\n" % (ii + 1, ' '.join([str(ii) for ii in x0]))
            # regions
            s += "# regions\n%d\n" % len(regions)
            for ii, x0 in enumerate(regions):
                s += " %d %f %f %d\n" % tuple([ii + 1] + x0)

        if self.dim == 3:

            for x in self.d2.values():
                assert isinstance(x, surface)
                p = [map[y.getn()] for y in x.getpoints()]
                h = []
                pts = []
                for hole in x.getholepoints():
                    h.append([map[y.getn()] for y in hole])
                    pts.append(getinsidepoint(hole).getxyz())
                bc = self.getBCnum(x.getn())
                facets.append((p, bc, h, pts))
            # number of facets, boundary markers=yes
            s += "# segments\n%d 1\n" % len(facets)
            for p, bc, h, holes in facets:
                # number of polygons, # of holes, boundary marker
                s += " %d %d %d\n" % (1 + len(h), len(h), bc)
                # number of corners, corner 1, corner 2, ...
                s += " %d %s\n" % (len(p), ' '.join([str(ii) for ii in p]))
                for x in h:
                    # number of corners, corner 1, corner 2, ...
                    s += " %d %s\n" % (len(x), ' '.join([str(ii) for ii in p]))
                for i, pt in enumerate(holes):
                    # hole #, x, y, z
                    s += ptstr % tuple([i + 1] + list(pt))

            # volume holes
            s += "# holes\n0\n"
            # regions
            regions=[]
            for x in self.phys3.values():
                assert isinstance(x, physicalvolume)
                for v in x.getvolumes():
                    regions.append(v.getinsidepoint().getxyz()+[x.getn()])
            s += "# regions\n%d\n" % len(regions)
            for i, x in enumerate(regions):
                s += ptstr2 % tuple([i + 1] + list(x))

        open(filename, "w").write(s)

    @staticmethod
    def from_gmsh_file(filename):
        """
        Import geometry - Gmsh geometry format.

        Parameters
        ----------
        filename : string
            file name

        Returns
        -------
        geo : geometry
            geometry description

        """

        from pyparsing import Word, Optional, nums, Combine, Literal, \
             CaselessLiteral, Group, OneOrMore, StringEnd, restOfLine, \
             ParseException, alphanums, Keyword, ZeroOrMore

        e = CaselessLiteral("E")
        inum = Word("+-"+nums)
        fnum = Combine(
            Word( "+-"+nums, nums ) + Optional("."+Optional(Word(nums))) +
            Optional(e+Word("+-"+nums,nums))
            )

        semi  = Literal(";").suppress()
        colon  = Literal(",").suppress()
        lpar  = Literal("(").suppress()
        rpar  = Literal(")").suppress()
        lbrace  = Literal("{").suppress()
        rbrace  = Literal("}").suppress()
        eq  = Literal("=").suppress()

        point = Group(
                Keyword("Point")+lpar+inum+rpar+eq+
                Group(lbrace+fnum+colon+fnum+colon+fnum+colon+fnum+rbrace)+semi
                )
        line = Group(
                Keyword("Line")+lpar+inum+rpar+eq+
                Group(lbrace+inum+colon+inum+rbrace)+semi
                )
        lineloop = Group(
                Keyword("Line Loop")+lpar+inum+rpar+eq+
                Group(lbrace+inum+OneOrMore(colon+inum)+rbrace)+semi
                )
        circle = Group(
                Keyword("Circle")+lpar+inum+rpar+eq+
                Group(lbrace+inum+colon+inum+colon+inum+rbrace)+semi
                )
        planesurface = Group(
                Keyword("Plane Surface")+lpar+inum+rpar+eq+
                Group(lbrace+inum+rbrace)+semi
                )
        ruledsurface = Group(
                Keyword("Ruled Surface")+lpar+inum+rpar+eq+
                Group(lbrace+inum+rbrace)+semi
                )
        surfaceloop = Group(
                Keyword("Surface Loop")+lpar+inum+rpar+eq+
                Group(lbrace+inum+OneOrMore(colon+inum)+rbrace)+semi
                )
        volume = Group(
                Keyword("Volume")+lpar+inum+rpar+eq+
                Group(lbrace+inum+rbrace)+semi
                )
        physicalsurface = Group(
                Keyword("Physical Surface")+lpar+inum+rpar+eq+
                Group(lbrace+inum+ZeroOrMore(colon+inum)+rbrace)+semi
                )
        physicalvolume = Group(
                Keyword("Physical Volume")+lpar+inum+rpar+eq+
                Group(lbrace+inum+ZeroOrMore(colon+inum)+rbrace)+semi
                )
        skip1 = Group(
                Word(alphanums)+eq+fnum+semi
                )

        comment = Group( Literal("//")+restOfLine).suppress()

        command = point | line | lineloop | circle | planesurface | ruledsurface | \
                surfaceloop | volume | physicalsurface | physicalvolume | comment \
                | skip1

        grammar= OneOrMore(command)+StringEnd()

        try:
            tokens= grammar.parseFile(filename)
        except ParseException as err:
            print(err.line)
            print(" "*(err.column-1) + "^")
            print(err)
            raise err

        lineloops={}
        surfaceloops={}
        geo=geometry()
        for x in tokens:
            if x[0]=="Point":
                geo.addpoint(int(x[1]),[float(x[2][0]),float(x[2][1]),float(x[2][2])])
            elif x[0]=="Line":
                assert len(x[2])==2
                geo.addline(int(x[1]),[int(x[2][0]),int(x[2][1])])
            elif x[0]=="Circle":
                assert len(x[2])==3
                geo.addline(int(x[1]),[int(x[2][0]),int(x[2][2])])
                #geo.add1(geom.circle(int(x[1]),int(x[2][0]),int(x[2][1]),
                #    int(x[2][2])))
            elif x[0]=="Line Loop":
                lineloops[int(x[1])]=[int(y) for y in x[2]]
            elif x[0]=="Plane Surface":
                assert len(x[2])==1
                geo.addsurface(int(x[1]),lineloops[int(x[2][0])])
            elif x[0]=="Ruled Surface":
                assert len(x[2])==1
                geo.addsurface(int(x[1]),lineloops[int(x[2][0])])
            elif x[0]=="Surface Loop":
                surfaceloops[int(x[1])]=[int(y) for y in x[2]]
            elif x[0]=="Volume":
                assert len(x[2])==1
                geo.addvolume(int(x[1]),surfaceloops[int(x[2][0])])
            elif x[0]=="Physical Surface":
                geo.addphysicalsurface(int(x[1]),[int(y) for y in x[2]])
            elif x[0]=="Physical Volume":
                geo.addphysicalvolume(int(x[1]),[int(y) for y in x[2]])
            else:
                raise "Unsupported entity: "+x[0]

        return geo

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
    def __truediv__(self,num):
        return self.__div__(num)
    def __mul__(self,num):
        return point(self.geom,-1,[a*num for a in self.p])
    def __rmul__(self,num):
        return self.__mul__(num)
    def getxyz(self):
        return self.p
    def getstr(self):
        if self.geom.dim == 2:
            return "%f, %f" % tuple(self.getxyz())
        elif self.geom.dim == 3:
            return "%f, %f, %f" % tuple(self.getxyz())
        else:
            return None

class line(geomobject):
    def __init__(self,g,n,l):
        self.geom=g
        self.n=n
        self.points=l
    def getpoints(self):
        return [self.geom.d0[x] for x in self.points]

class surface(geomobject):
    def __init__(self,g,n,s, is_hole=False):
        self.geom=g
        self.n=n
        self.lines,self.holes=self.separate(s)
        self.is_hole = is_hole
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
        if self.is_hole:
            r.append(self.getinsidepoint())
        else:
            for hole in self.holes:
                r.append([p(x) for x in hole])
        return r
    def getcenterpoint(self):
        pts=self.getpoints()
        return sum(pts[1:], pts[0] * 0) / float(len(pts))
    def getinsidepoint(self):
        p0 = self.getpoints()[0]
        pc = self.getcenterpoint()
        return p0 + (pc - p0) * 0.001

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
        direct=(pts[0]+pts[1]+pts[2]) / 3.0 - p0
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

import geometry as geom

def getinsidepoint(pts):
    direct = (pts[0] + pts[1] + pts[2]) / 3 - pts[0]
    return pts[0] + 0.001 * direct

def write_poly(g, filename):
    if g.dim == 2:
        g.leaveonlyphysicalsurfaces()
    if g.dim == 3:
        g.leaveonlyphysicalvolumes()

    # write nodes
    nodes = []
    map = {}
    for x in g.d0.values():
        assert isinstance(x, geom.point)
        nodes.append(x.getxyz())
        map[x.getn()] = len(nodes)

    s = "# nodes\n%d %d 0 0\n" % (len(nodes), g.dim)
    if g.dim == 2:
        ptstr = " %d %f %f\n"
    else:
        ptstr = " %d %f %f %f\n"

    for n, x in enumerate(nodes):
        s += ptstr % tuple([n + 1] + list(x))

    # facets
    # first write external polygon, then hole polygons and then point in each
    # hole polygon
    facets = []
    if g.dim == 2:

        hole_pts = []
        regions=[]
        for x2 in g.d2.values():
            assert isinstance(x2, geom.surface)
            for x1 in x2.getlines():
                assert isinstance(x1, geom.line)
                p = [map[y.getn()] for y in x1.getpoints()]
                bc = g.getBCnum(x1.getn())
                facets.append((p, bc))

            for hole in x2.getholepoints():
                hole_pts.append(hole.getxyz())

        # regions
        for x in g.phys2.values():
            assert isinstance(x, geom.physicalsurface)
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

    if g.dim == 3:

        for x in g.d2.values():
            assert isinstance(x, geom.surface)
            p = [map[y.getn()] for y in x.getpoints()]
            h = []
            pts = []
            for hole in x.getholepoints():
                h.append([map[y.getn()] for y in hole])
                pts.append(getinsidepoint(hole).getxyz())
            bc = g.getBCnum(x.getn())
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
        for x in g.phys3.values():
            assert isinstance(x, geom.physicalvolume)
            for v in x.getvolumes():
                regions.append(v.getinsidepoint().getxyz()+[x.getn()])
        s += "# regions\n%d\n" % len(regions)
        for i, x in enumerate(regions):
            s += ptstr % tuple([i + 1], list(x))

    open(filename, "w").write(s)

def read_poly_mesh(fname, verbose=True):

    def getnodes(fnods,up):
        f=file(fnods)
        l=[int(x) for x in f.readline().split()]
        npoints,dim,nattrib,nbound=l
        if verbose: up.init(npoints)
        nodes=[]
        for line in f:
            if line[0]=="#": continue
            l=[float(x) for x in line.split()]
            l = l[:(dim + 1)]
            l[0]=int(l[0])
            nodes.append(tuple(l))
            assert l[0]==len(nodes)
        assert npoints==len(nodes)
        return nodes

    def getele(fele,up):
        f=file(fele)
        l=[int(x) for x in f.readline().split()]
        nele,nnod,nattrib=l
        #we have either linear or quadratic tetrahedra:
        if nnod in [4,10]:
            elem = 'tetra'
            linear = (nnod == 4)
        if nnod in [3, 7]:
            elem = 'tri'
            linear = (nnod == 3)

        # if nattrib!=1:
        #     raise "tetgen didn't assign an entity number to each element (option -A)"
        els=[]
        regions={}
        for line in f:
            if line[0]=="#": continue
            l=[int(x) for x in line.split()]
            if elem == 'tri':
                if linear:
                    assert (len(l) - 1 - nattrib) == 3
                    els.append((l[0],l[1],l[2],l[3]))
                    regionnum=l[5]
                else:
                    assert len(l)-2 == 10
                    els.append((l[0],54,l[1],l[2],l[3],l[4],
                                l[5],l[6],l[7],l[8],l[9],l[10]))
                    regionnum=l[11]
            if elem == 'tetra':
                if linear:
                    assert len(l)-2 == 4
                    els.append((l[0],54,l[1],l[2],l[3],l[4]))
                    regionnum=l[5]
                else:
                    assert len(l)-2 == 10
                    els.append((l[0],54,l[1],l[2],l[3],l[4],
                                l[5],l[6],l[7],l[8],l[9],l[10]))
                    regionnum=l[11]
            if regionnum==0:
                print "see %s, element # %d"%(fele,l[0])
                raise "there are elements not belonging to any physical entity"
            if regions.has_key(regionnum):
                regions[regionnum].append(l[0])
            else:
                regions[regionnum]=[l[0]]
            assert l[0]==len(els)
            if verbose: up.update(l[0])
        return els,regions,linear

    def getBCfaces(ffaces,up):
        f=file(ffaces)
        l=[int(x) for x in f.readline().split()]
        nfaces,nattrib=l
        if nattrib!=1:
            raise "tetgen didn't assign an entity number to each face \
(option -A)"
        if verbose: up.init(nfaces)
        faces={}
        for line in f:
            if line[0]=="#": continue
            l=[int(x) for x in line.split()]
            assert len(l)==5
            regionnum=l[4]
            if regionnum==0: continue
            if faces.has_key(regionnum):
                faces[regionnum].append((l[1],l[2],l[3]))
            else:
                faces[regionnum]=[(l[1],l[2],l[3])]
            if verbose: up.update(l[0])
        return faces

    def calculatexyz(nodes, els):
        """Calculate the missing xyz values in place"""
        def avg(i,j,n4,nodes):
            a=nodes[n4[i-1]-1]
            b=nodes[n4[j-1]-1]
            return (a[1]+b[1])/2, (a[2]+b[2])/2, (a[3]+b[3])/2
        def getxyz(i,n4,nodes):
            if i+5==5: return avg(1,2,n4,nodes)
            if i+5==6: return avg(2,3,n4,nodes)
            if i+5==7: return avg(1,3,n4,nodes)
            if i+5==8: return avg(1,4,n4,nodes)
            if i+5==9: return avg(2,4,n4,nodes)
            if i+5==10: return avg(3,4,n4,nodes)
            raise "wrong topology"
        for e in els:
            n4=e[2:2+4]
            n6=e[2+4:2+4+10]
            for i,n in enumerate(n6):
                x,y,z=getxyz(i,n4,nodes)
                nodes[n-1]=(n,x,y,z)

    if verbose: print "Reading geometry from poly file..."
    m=mesh()
    m.nodes=getnodes(fname+".node")
    m.elements,m.regions, lin=getele(fname+".ele")
    if not lin:
        #tetgen doesn't compute xyz coordinates of the aditional 6 nodes
        #(only of the 4 corner nodes) in tetrahedra.
        calculatexyz(m.nodes,m.elements)
    m.faces=getBCfaces(fname+".face")
    return m

def read_gmsh(filename):
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
    except ParseException, err:
        print err.line
        print " "*(err.column-1) + "^"
        print err
        raise err

    lineloops={}
    surfaceloops={}
    g=geom.geometry()
    for x in tokens:
        if x[0]=="Point":
            g.addpoint(int(x[1]),[float(x[2][0]),float(x[2][1]),float(x[2][2])])
        elif x[0]=="Line":
            assert len(x[2])==2
            g.addline(int(x[1]),[int(x[2][0]),int(x[2][1])])
        elif x[0]=="Circle":
            assert len(x[2])==3
            g.addline(int(x[1]),[int(x[2][0]),int(x[2][2])])
            #g.add1(geom.circle(int(x[1]),int(x[2][0]),int(x[2][1]),
            #    int(x[2][2])))
        elif x[0]=="Line Loop":
            lineloops[int(x[1])]=[int(y) for y in x[2]]
        elif x[0]=="Plane Surface":
            assert len(x[2])==1
            g.addsurface(int(x[1]),lineloops[int(x[2][0])])
        elif x[0]=="Ruled Surface":
            assert len(x[2])==1
            g.addsurface(int(x[1]),lineloops[int(x[2][0])])
        elif x[0]=="Surface Loop":
            surfaceloops[int(x[1])]=[int(y) for y in x[2]]
        elif x[0]=="Volume":
            assert len(x[2])==1
            g.addvolume(int(x[1]),surfaceloops[int(x[2][0])])
        elif x[0]=="Physical Surface":
            g.addphysicalsurface(int(x[1]),[int(y) for y in x[2]])
        elif x[0]=="Physical Volume":
            g.addphysicalvolume(int(x[1]),[int(y) for y in x[2]])
        else:
            raise "Unsupported entity: "+x[0]
    return g

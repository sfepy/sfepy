import math

import geometry as geom

from meshutils import mesh
import progressbar

def numlist2str(x):
    s=""
    for i in x:
        s+="%d "%i
    return s[:-1]

def getinsidepoint(pts):
    direct=(pts[0]+pts[1]+pts[2])/3-pts[0]
    return pts[0]+0.001*direct

def write_tetgen(g,filename):
    g.leaveonlyphysicalvolumes()
    #nodes
    nodes=[]
    map={}
    for x in g.d0.values():
        assert isinstance(x,geom.point)
        nodes.append(x.getxyz())
        map[x.getn()]=len(nodes)
    s="%d 3\n"%len(nodes)
    for n,x in enumerate(nodes):
        s+="%d %f %f %f\n"%tuple([n+1]+list(x))

    #facets
    #first write external polygon, then hole polygons and then point in each
    #hole polygon
    facets=[]
    for x in g.d2.values():
        assert isinstance(x,geom.surface)
        p=[map[y.getn()] for y in x.getpoints()]
        h=[]
        pts=[]
        for hole in x.getholepoints():
            h.append([map[y.getn()] for y in hole])
            pts.append(getinsidepoint(hole).getxyz())
        bc=g.getBCnum(x.getn())
        facets.append((p,bc,h,pts))
    # # of facets, boundary markers=yes
    s+="\n%d 1\n"%len(facets)
    for p,bc,h,holes in facets:
        # # of polygons, # of holes, boundary marker
        s+="%d %d %d\n"%(1+len(h),len(h),bc)
        # # of corners, corner 1, corner 2, ...
        s+="%d %s\n"%(len(p),numlist2str(p))
        for x in h:
            # # of corners, corner 1, corner 2, ...
            s+="%d %s\n"%(len(x),numlist2str(x))
        for i,pt in enumerate(holes):
            # hole #, x, y, z
            s+="%d %f %f %f\n"%(i+1,pt[0],pt[1],pt[2])

    #volume holes
    s+="\n0\n"

    #regions
    regions=[]
    for x in g.phys3.values():
        assert isinstance(x,geom.physicalvolume)
        for v in x.getvolumes():
            regions.append(v.getinsidepoint().getxyz()+[x.getn()])
    s+="\n%d\n"%len(regions)
    for i,x in enumerate(regions):
        s+="%d %f %f %f %d\n"%(i+1,x[0],x[1],x[2],x[3])
    open(filename,"w").write(s)

def read_tetgen(fname,verbose=True):
    def getnodes(fnods,up):
        f=file(fnods)
        l=[int(x) for x in f.readline().split()]
        npoints,dim,nattrib,nbound=l
        assert dim==3
        if verbose: up.init(npoints)
        nodes=[]
        for line in f:
            if line[0]=="#": continue
            l=[float(x) for x in line.split()]
            l[0]=int(l[0])
            nodes.append(tuple(l))
            assert l[0]==len(nodes)
            if verbose: up.update(l[0])
        assert npoints==len(nodes)
        return nodes
    def getele(fele,up):
        f=file(fele)
        l=[int(x) for x in f.readline().split()]
        ntetra,nnod,nattrib=l
        #we have either linear or quadratic tetrahedra:
        assert nnod in [4,10]
        linear= (nnod==4)
        if verbose: up.init(ntetra)
        if nattrib!=1:
            raise "tetgen didn't assign an entity number to each element \
(option -A)"
        els=[]
        regions={}
        for line in f:
            if line[0]=="#": continue
            l=[int(x) for x in line.split()]
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

    if verbose: print "Reading mesh from tetgen..."
    m=mesh()
    m.nodes=getnodes(fname+".node",progressbar.MyBar("        nodes:"))
    m.elements,m.regions, lin=getele(fname+".ele",
            progressbar.MyBar("        elements:"))
    if not lin:
        #tetgen doesn't compute xyz coordinates of the aditional 6 nodes
        #(only of the 4 corner nodes) in tetrahedra.
        calculatexyz(m.nodes,m.elements)
    m.faces=getBCfaces(fname+".face",progressbar.MyBar("        BC:"))
    return m

def runtetgen(filename,a=None,Q=None,quadratic=False,verbose=True,
        refine=False,tetgenpath="/usr/bin/tetgen"):
    """Runs tetgen.
    
    tetgenpath ... the tetgen executable with a full path
    filename ... the input file for tetgen (for example /tmp/t.poly)
    a ... a maximum tetrahedron volume constraint
    Q ... a minimum radius-edge ratio, tetgen default is 2.0
    quadratic ... False - generate linear elements, True - quadratic elements
    """
    import pexpect
    if not refine:
        cmd = "%s -pQAq" % (tetgenpath)
    else:
        cmd = "%s -rQAq" % (tetgenpath)
    if Q!=None:
        cmd=cmd+"%f"%Q
    if a!=None and not refine:
        cmd=cmd+" -a%f"%(a)
    if refine:
        cmd=cmd+" -a"
    if quadratic:
        cmd=cmd+" -o2"
    cmd=cmd+" %s"%(filename)
    if verbose: print "Generating mesh using", cmd
    p=pexpect.spawn(cmd,timeout=None)
    if not refine:
        p.expect("Opening %s."%(filename))
    else:
        p.expect("Opening %s.node.\r\n"%(filename))
        p.expect("Opening %s.ele.\r\n"%(filename))
        p.expect("Opening %s.face.\r\n"%(filename))
        p.expect("Opening %s.vol."%(filename))
    assert p.before==""
    p.expect(pexpect.EOF)
    if p.before!="\r\n":
        print p.before
        raise "Error when running tetgen (see above for output): %s"%cmd

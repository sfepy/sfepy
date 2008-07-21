import math

import geometry as geom

def conv(x):
    s="["
    for i in x:
        s+=str(i)+","
    return s[:-1]+"]"

def conv2(x):
    s=""
    for i in x:
        s+=i+","
    return s[:-1]

def conv3(x):
    s=""
    for i in x:
        s+="'"+i+"',"
    return s[:-1]

def conv4(a,b,c):
    s="["
    for i in a:
        s+=str(i)+" "
    s=s[:-1]+"; "
    for i in b:
        s+=str(i)+" "
    s=s[:-1]+"; "
    for i in c:
        s+=str(i)+" "
    return s[:-1]+"]"


def axissym(a,b,c):
    "a,b ... line. c ... center"

    d=[(ai+bi)/2.-ci for ai,bi,ci in zip(a,b,c)]
    dn=math.sqrt(d[0]**2+d[1]**2+d[2]**2)
    return [ci+di+0.1*di/dn for ci,di in zip(c,d)]

def femlabsurface3_old(f,n,p):
    #fucking femlab - there is a bug in face3 for some points
    assert len(p)==3
    a=[y[0] for y in p]
    b=[y[1] for y in p]
    c=[y[2] for y in p]
    return (f+"=face3(%s',%s',%s')\n")%(n,conv(a),conv(b),conv(c))

def curve2(f,p1,p2):
    return "%s=curve2(%s,%s);\n"%(f,conv([p1[0],p2[0]]),conv([p1[1],p2[1]]))

def getp3(p1,p2,p3):
    import math
    a=norm2(vec(p2,p3))
    b=norm2(vec(p1,p3))
    c=norm2(vec(p1,p2))
    cosphi=(b**2+c**2-a**2)/(2*b*c)
    if cosphi>1.0: cosphi=1.0
    return [b*cosphi, b*math.sqrt(1-cosphi**2)]


def femlabsurface3(f,n,p):
    assert len(p)==3
    variable=f%n
    f1=variable+"f1"
    f2=variable+"f2"
    f3=variable+"f3"
    g4=variable+"g4"
    s=""
    p1=[0,0]
    p2=[norm2(vec(p[1],p[0])),0]
    p3=getp3(p[0],p[1],p[2])
    s+=curve2(f1,p1,p2)
    s+=curve2(f2,p2,p3)
    s+=curve2(f3,p3,p1)
    s+="%s=geomcoerce('solid',{%s,%s,%s});\n"%(g4,f1,f2,f3)
    a=[y[0] for y in p]
    b=[y[1] for y in p]
    c=[y[2] for y in p]
#    print p,norm2(crossproduct(vec(p[0],p[1]), vec(p[0],p[2])))
    s+="%s=embed(%s,'Wrkpln',%s);\n"%(variable,g4,conv4(a,b,c))
    return s

def femlabline(f,n,p):
    assert len(p)==2
    a=[y.getxyz()[0] for y in p]
    b=[y.getxyz()[1] for y in p]
    c=[y.getxyz()[2] for y in p]
    return (f+"=curve3(%s,%s,%s)\n")%(n,conv(a),conv(b),conv(c))

def triangulate2(points):
    if len(points)==3:
        return (points,)
    if len(points)==4:
        return ( (points[0],points[1],points[2]),
                (points[0],points[2],points[3]) )
    import Polygon 
    q=Polygon.Polygon(points)
    strips=Polygon.TriStrip(q)
    tri=[]
    for t in strips:
        for n in range(len(t)-2):
            tri.append((t[n],t[n+1],t[n+2]))
    return tri

def triangulate(points):
    from poly import poly
#    if len(points)>3:
#        print "-"*60
#        print points
#        print poly.triangulate(points)
#        print "-"*60
    return poly.triangulate(points)

def crossproduct(a,b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

def norm2(x):
    import math
    return math.sqrt(x[0]**2+x[1]**2+x[2]**2)

def norm(x):
    n=norm2(x)
    return [xi/n for xi in x]

def vec(a,b):
    return [ai-bi for ai,bi in zip(a,b)]

def getnormal(a,b,c):
    v1=[ai-bi for ai,bi in zip(a,b)]
    v2=[ci-ai for ai,ci in zip(a,c)]
    return norm(crossproduct(v1,v2))

def tri3d(points):
    n=getnormal(points[0],points[1],points[2])
    p=points[0]
    d=-(n[0]*p[0]+n[1]*p[1]+n[2]*p[2])
    #print n,d
    for p in points:
        f=n[0]*p[0]+n[1]*p[1]+n[2]*p[2]+d
        assert abs(f)<1e-8
    z=0
    if abs(n[1])>abs(n[z]):z=1
    if abs(n[2])>abs(n[z]):z=2
    x=0
    if x==z: x=1
    y=1
    if y==x or y==z: y=2
    p2=[(p[x],p[y]) for p in points]
    t=triangulate(p2)
    result=[]
    for i in range(len(t)):
        pp=[]
        for j in range(3):
            p=[0,0,0]
            p[x]=t[i][j][0]
            p[y]=t[i][j][1]
            p[z]=-(n[x]*p[x]+n[y]*p[y]+d)/n[z]
            pp.append(p)
        result.append(pp)
    return result

def triangulate_old(p):
    tri=[]
    for i in range(1,len(p)-1):
        tri.append((p[0],p[i],p[i+1]))
    return tri

def femlabsurface(f,n,points):
    s=""
    x=f+"=geomcoerce('face',{"
    p=[y.getxyz() for y in points]
#    print n,":",len(p)
    tri=tri3d(p)
    for i,t in enumerate(tri):
        if len(tri)==1:
            ff=f
        else:
            ff=f+"P%d"%i
            x+=f+"P%d,"%i
        s+=femlabsurface3(ff,n,t)
    if len(tri)>1:
        x=x[:-1]+"})\n"
        s+=x%tuple([n]*(len(tri)+1))
    return s


def write_femlab(g,filename, export0D=False, export1D=False, export2D=False,
    export3D=False):
    if not export1D and not export2D and not export3D and not export3D:
        export3D=True
    head="""\
flclear fem

% COMSOL version
clear vrsn
vrsn.name = 'COMSOL 3.2';
vrsn.ext = '';
vrsn.major = 0;
vrsn.build = 222;
vrsn.rcs = '$Name:  $';
vrsn.date = '$Date: 2005/09/01 18:02:30 $';
fem.version = vrsn;

"""
    tail="""
clear appl
appl.mode.class = 'ConductiveMediaDC';
appl.shape = {};
appl.gporder = {};
appl.cporder = {};
appl.sshape = 2;
appl.assignsuffix = '_dc';
clear pnt
pnt.V0 = {};
pnt.type = {};
pnt.Qj0 = {};
pnt.name = {};
pnt.ind = [];
appl.pnt = pnt;
clear edg
edg.Qlj = {};
edg.name = {};
edg.ind = [];
appl.edg = edg;
clear bnd
bnd.Vref = {};
bnd.sigmabnd = {};
bnd.V0 = {};
bnd.Jn = {};
bnd.type = {};
bnd.dbnd = {};
bnd.J0 = {};
bnd.name = {};
bnd.ind = [];
appl.bnd = bnd;
clear equ
equ.init = {};
equ.cporder = {};
equ.T0 = {};
equ.res0 = {};
equ.gporder = {};
equ.Qj = {};
equ.sigma = {};
equ.usage = {};
equ.T = {};
equ.name = {};
equ.Je = {};
equ.sigmatensor = {};
equ.sigtype = {};
equ.alpha = {};
equ.ind = [];
appl.equ = equ;

fem.appl{1} = appl;
fem.sdim = {'x','y','z'};
fem.border = 1;
fem.units = 'SI';

% Multiphysics
fem=multiphysics(fem);
"""
    s=""
    objs=[]
    if export0D:
        for x in g.d0.values():
            assert isinstance(x,geom.point)
            s+="p%d=point3(%s)\n"%(x.getn(),x.getstr())
            objs.append("p%d"%x.getn())
    if export1D:
        for x in g.d1.values():
            if isinstance(x,geom.line):
                p=x.getpoints()
                s+=femlabline("l%s",x.getn(),p)
            elif isinstance(x,geom.circle):
                p=x.getpoints()
                assert len(p)==3
                s+=femlabline("l%s",x.getn(),(p[0],p[2]))
            elif isinstance(x,geom.lineloop):
                continue
            else:
                print "Warning: unknown element ",type(x)
                continue
            objs.append("l%d"%x.getn())
    if export2D or export3D:
        for x in g.d2.values():
            if isinstance(x,geom.planesurface):
                p=x.getpoints()
                s+=femlabsurface("f%s",x.getn(),p);
            elif isinstance(x,geom.ruledsurface):
                p=x.getpoints()
                s+=femlabsurface("f%s",x.getn(),p);
            elif isinstance(x,geom.surfaceloop):
                continue
            else:
                print "Warning: unknown element ",type(x)
                continue
            if export2D: 
                objs.append("f%s"%x.getn())
    if export3D:
        for x in g.d3.values():
            if isinstance(x,geom.volume):
                p=x.getsurfaces()
                s+="s%s"%x.getn()+"=geomcoerce('solid',{"
                for y in p:
                    s+="f%s,"%y.getn()
                s=s[:-1]+"})\n"
            else:
                print "Warning: unknown element ",type(x)
                continue
            objs.append("s%s"%x.getn())

    s+="clear s\ns.objs={%s};\ns.name={%s};\nfem.draw=struct('s',s);\n"%\
        (conv2(objs),conv3(objs))
    s= head+s+tail
    open(filename,"w").write(s)

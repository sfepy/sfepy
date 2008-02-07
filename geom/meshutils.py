# meshutils.py
"""
Finite element mesh utilites.
"""
__docformat__ = "restructuredtext en"

# Copyright (C) 2004-2005 O. Certik 
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA.
#
# Please report all bugs and problems to <ondrej@certik.cz>.

################################################

import string
import math
import os

from pyparsing import Word, Optional, alphas, nums, Combine, Literal, CaselessLiteral, LineEnd, Group, Dict, OneOrMore, StringEnd, restOfLine, ParseException, oneOf, Forward, alphanums

import progressbar

#gmsh element types, see 
#http://www.geuz.org/gmsh/doc/texinfo/gmsh_10.html#SEC65
#1D elements
mshpoint=15 #Point (1 node)
mshline=1 #Line (2 nodes)
mshline2=8 #Second order line (3 nodes)
#2D elements:
mshtriangle=2 #Triangle (3 nodes)
mshtriangle2=9 #Second order triangle (6 nodes)
mshquadrangle=3 #Quadrangle (4 nodes)
mshquadrangle2=10 #Second order quadrangle (9 nodes)
#3D elements
mshtetrahedron=4 #Tetrahedron (4 nodes)
mshtetrahedron2=11 #Tetrahedron (10 nodes)
mshprism=6 #Prism (6 nodes)
mshhexahedron=5 #Hexahedron (8 nodes)

#pmd element types, see
#http://www.it.cas.cz/manual/pmd/ra.htm
#2D elements:
pmdtriangle=4 #Triangle (3 or 6 nodes)
pmdtrianglerot=5 #Triangle (3 or 6 nodes)
pmdquadrangle=6 #Quadrangle (4 or 8 nodes)
pmdquadranglerot=7 #Quadrangle (4 or 8 nodes)
#3D elements
pmdtetrahedron=54 #Tetrahedron (4 or 10 nodes)
pmdprism=55 #Prism (6 or 15 nodes)
pmdhexahedron=56 #Hexahedron (8 or 20 nodes)
pmdsemiloof=61 #triangular semi-loof (6 nodes)

#libmesh element types, see
#http://libmesh.sourceforge.net/doxygen/namespacelibMeshEnums.php#54ee290fca7f0c26eb1e5986f52714574518904c8c2948ef3d2869c4cb4a2b8f
#2D elements:
libmeshtriangle=3 #TRI3
libmeshquadrangle=5 #QUAD4
#3D elements
libmeshtetrahedron=8 #TET4
libmeshtetrahedron_quad=9 #TET10
libmeshhexahedron=10 #HEX8
libmeshprism=13 #PRISM6

#number of the physical entity, which represents the whole model
mshmodelnum=100

def check(s,what):
    if s != what: 
        error("'%s' missing"%(what),1)

class MeshUtilsError(Exception):
    """Base class for exceptions in meshutils."""
    pass

class MeshUtilsCheckError(MeshUtilsError):
    """Function check failed."""
    pass

class MeshUtilsParseError(MeshUtilsError):
    """Parse error in input file."""
    pass
class MeshUtilsWarning(MeshUtilsError):
    """Not necessarily error, but some strange thing happened."""
    pass

def error(s,type=0):
    if type==1:
        raise MeshUtilsCheckError,s
    elif type==2:
        raise MeshUtilsParseError,s
    elif type==3:
        raise MeshUtilsWarning,s
    else:
        raise MeshUtilsError,s

def myfloat(s):
    "Converts s to float, including PMD float format (without E)."
    try:
        f=float(s)
    except ValueError:
        s=s[:-2]+"E"+s[-2:]
        f=float(s)
    return f
    


class bound:
    """Handles all physical entities <= mshmodelnum.

    This class is optionaly filled in mesh.readmsh() by calling the method
    bound.handle2(). It extracts **only** the node numbers of the entity.
    All entities are stored in self.f dictionary (key is the entity number).

    Has methods to get the node numbers in PMD notation:
    (1 2 4 1 2 3 4 8 9 10 -> 1:4 8:10)
    ie. it removes any repeating numbers and shortens the list using ":" 

    Has method to find elements from a given list, which lie on the
    boundary formed by nodes (also returns appropriate side of elements).

    So - all boundary conditions should be done using this class.
    """
    def __init__(self):
        self.f={} #dictionary for node numbers
        self.elementsassociated=False
    def readpmd(self,filename):
        """Loads internal dictionary from 'filename'.

        Format must be the same as from write()"""
        self.f={}
        f=file(filename,"r")
        line=f.readline()
        while line:
            x=string.split(line)
            n=int(x[1])
            if self.f.has_key(n):
                error("two /N have the same number!",3)
            self.f[n]=self.fromstr(x[3:len(x)])
            line=f.readline()
        f.close()

    def fromstr(self,str):
        f=[]
        try:
            for x in str:
                seq=string.split(x,":")
                if len(seq) == 1:
                    f.append(int(x))
                elif len(seq) == 2:
                    f.extend(range(int(seq[0]),
                        int(seq[1])+1))
                else:
                    error("Invalid syntax",2)
        except ValueError:
            error("Invalid syntax",2)
        return f
    def writepmd(self,filename):
        f=file(filename,"w")
        for n,l in self.f.iteritems():
            f.write("  /N %d N"%(n))
            f.write(self.str(self.simplify(l)))
            f.write("\n")
        f.close()
    def getstr(self,key):
        return self.str(self.getf(key))
    def getf(self,key):
        if key > mshmodelnum and not self.elementsassociated:
            error("Elements from entity %d aren't associated."%(key))
        return self.simplify(self.f[key])
    def simplify(self,l):
        """Removes repeating numbers and sorts the internal list l.
        
        Note: it's very slow for large meshes. This is the thing
        which slows down everything.
        """
        q=[]
        for x in l:
            if not (x in q): q.append(x)
        q.sort()
        return q
    def str(self,l):
        "Converts l (must be a sorted list) to string (using ':')."
        s=""
        old2x=-1
        oldx=-1
        for x in l: 
            if x!=oldx+1:
                if (oldx != -1) and (oldx!=old2x):
                    if oldx==old2x+1:
                        s+=" %d"%(oldx)
                    else:
                        s+=":%d"%(oldx)
                s+=" %d"%(x)
                old2x=x
            oldx=x
        if (oldx != -1) and (oldx!=old2x):
            if oldx==old2x+1:
                s+=" %d"%(oldx)
            else:
                s+=":%d"%(oldx)

        return s[1:]
    def handle(self,n,list):
        "Appends number in 'list' to internal dictionary (key=n)"
        if not self.f.has_key(n):
            self.f[n]=[]
        self.f[n].extend(list)
    def handle2(self,p):
        entity=p[2]
        eltype=p[1]
        if eltype==mshpoint:
            self.handle(entity,(p[5],))
        elif eltype==mshline:
            self.handle(entity,(p[5],p[6]))
        elif eltype==mshline2:
            self.handle(entity,(p[5],p[6],p[7]))
        elif eltype==mshtriangle:
            self.handle(entity,(p[5],p[6],p[7]))
        elif eltype==mshquadrangle:
            self.handle(entity,(p[5],p[6],p[7],p[8]))
        else:
            error("unsupported element type. "\
                "entity: %d; eltype %d;\n%s"\
                %(entity,eltype,repr(p)),3)
    def handleelement(self,p):
        """Handles element - not for boundary conditions"""
        self.elementsassociated=False
        entity=p[2]
        if not self.f.has_key(entity):
            self.f[entity]=[]
        self.f[entity].append([p[0],p[1]]+p[5:])
    def findelements(self,key,elements):
        """Returns element numbers (in a list), which lies
        at the boundary determined by self.nodes.
        
        Also returns the side of element. Currently only support
        triangles and quadrangles.
        """
        if not self.f.has_key(key):
            error("physical entity %d isn't in bound, aborting..."%(
                key))
        el=[]
        nodes=self.simplify(self.f[key])
        #print nodes
        for p in elements:
            nods=[]
            i=1
            for n in p[2:]: 
                if nodes.count(n): 
                    nods.append(i)
                i+=1
            if self.is2d:
                if len(nods)==3: 
                    if p[1] == pmdtriangle or \
                        p[1] == pmdtrianglerot:
                        while nods[-1] > 3:
                            nods=nods[:-1]
                    if p[1] == pmdquadrangle or \
                        p[1] == pmdquadrangle:
                        while nods[-1] > 4:
                            nods=nods[:-1]
                if len(nods)==2: 
                    if nods==[1,2]:
                        side=1
                    elif nods==[2,3]:
                        side=2 
                    elif nods==[3,4]:
                        side=3 
                    elif nods==[1,3]:
                        side=3 #must be triangle. check it.
                        if p[1] != pmdtriangle and p[1] != pmdtrianglerot:
                            error("findelements: error 3",3)
                    elif nods==[1,4]:
                        side=4 #must be quadrangle. check it.
                        if p[1] != pmdquadrangle and p[1] != pmdquadranglerot:
                            error("findelements: error 4",3)
                    else:
                        error("findelements: error 2",3)
                    el.append((p[0],side))
            else:
                if len(nods) == 3:
                    if p[1] == pmdtetrahedron:
                        if nods==[1,2,3]:
                            side=1
                        elif nods==[1,2,4]:
                            side=2
                        elif nods==[2,3,4]:
                            side=3
                        elif nods==[1,3,4]:
                            side=4
                        else:
                            side=0
                            error("findelements:error 7",3)
                        el.append((p[0],side))
                    else:
                        error("findelements:error 6",3)
                if len(nods)==4:
                    if p[1] == pmdhexahedron:
                        if nods==[1,2,3,4]:
                            side=1
                        elif nods==[1,2,5,6]:
                            side=2
                        elif nods==[2,3,6,7]:
                            side=3
                        elif nods==[3,4,7,8]:
                            side=4
                        elif nods==[1,4,5,8]:
                            side=5
                        elif nods==[5,6,7,8]:
                            side=6
                        else:
                            side=0
                            error("findelements:error 7",3)
                        el.append((p[0],side))
                    elif p[1] == pmdprism:
                        if nods == [1,3,4,6]:
                            side=4
                        else:
                            side=0
                            error("findelements:error 9",3)
                        el.append((p[0],side))
                    else:
                        error("findelements:error 5",3)
        return el
    def writeSV(self,f,els,key):
        for e in els:
            f.write("  /S %d E %d S%d\n"%(key,e[0],e[1]))
    def associateelements(self,elements):
        if self.elementsassociated:
            return
        keys=[k for k in self.f.keys() if k>mshmodelnum]
        for key in keys:
            self.associateelements_key_fast(key,elements)
        self.elementsassociated=True
    def associateelements_key_fast(self,key,elements):
        """Finds the elements under the key ``key`` in ``elements``.

        This is a fast version, which assumes an undocumented feature, that
        all the elements which gmsh exports are in *exactly* the same
        order both in the entity `key` and `mshmodelnum`.
        """
        list=self.f[key]
        p=list[0]
        first=self.finde(p,elements)
        listout=range(first,first+len(list))
        #print listout
        #listout.sort()
        #print listout
        self.f[key]=listout
    def associateelements_key(self,key,elements):
        """Finds the elements under the key ``key`` in ``elements``.

        This is a regular (very slow, but bullet proof) version.
        """
        list=self.f[key]
        listout=[]
        for n,p in enumerate(list):
            listout.append(self.finde(p,elements))
        #listout.sort()
        #print listout
        self.f[key]=listout
    def finde(self,p,elements):
        elnum=0
        for e in elements:
            if tuple(e[2:])==tuple(p[2:]):
                elnum=e[0]
                break
        if elnum==0:
            error("element not found.",2)
        return elnum

class mesh:
    #This class should be refactored anyway, to include better IO
    #interface, as I was thinking a year ago
    #the bound() is probably unnecessary, and the IO routines 
    #(pmd,libmesh,gmsh,tetgen...) should be done in some clever way
    #currently tetgen is in the module tetgen, others are here, then
    #the geometry is in tetgen and gmsh modules... 
    #ideal solution is to have 2 classes: geometry and mesh. both would
    #support import/export in some clean way.
    def __init__(self):
        self.clean()
    def clean(self):
        "Deletes the whole mesh."
        #following variables are the only variables in mesh 
        #(=the state mesh instance is in is fully determined be
        #these variables)
        #the user can do what he wants with them
        #and wherever the method mesh.clean() is called,
        #the mesh instance is "fresh" afterwards.
        #also it isn't allowed to store any other variable in
        #mesh class then the following.

        self.nodes = [] 
        #(n,x,y,z) : (int,float,float,float)
        #n ... node number
        #  always ordered in a consecutive way, starting from 1
        #  it must already be given like this in I1 and msh
        #x,y,z ... node coordinates
        self.elements = [] 
        #(n, eltype, nodes....) : (int,int,int,int,int...) 
        #n .... element number
        #  always ordered in a consecutive way, starting from 1
        #  it must already be given like this in I1 and msh
        #eltype .... PMD element type
        #nodes ... associated node numbers
        #   1st and 2nd order elements differs only by the number of
        #   nodes (ie. 3 and 6 for a triangle)
        #
        #The whole mesh is in fact PMD type of mesh - it is using
        #PMD primitives etc. So conversion to/from other formats
        #is handled by appropriate functions (writemsh,readmsh,
        #readELE,writeELE,...)
        self.boundbox = (0,0,0,0,0,0) 
        #reading:
        #  I1: read from RP
        #  msh: computed from node coordinates and 
        #     added +-eps (see readmsh)
        #  ELE,NOD:isn't set
        #writing:
        #  I1: written to RP
        #  msh,ELE,NOD: isn't used
        self.is2d=False #is the problem 2D (z==0) ?
        #reading:
        #  I1: set automatically distinguishing between (x,y) and(x,y,z)
        #  msh: set automatically dis. (x,y,"0") and (x,y,z)
        #  NOD: set automatically dis. (x,y,"0.0") and (x,y,z)
        #  ELE: isn't set
        #writing:
        #  I1:(x,y) versus (x,y,z)
        #  msh:(x,y,"0") versus (x,y,z)
        #  NOD,ELE: isn't used
        self.symmetric = False #rotational symmetric
        #reading:
        #  I1: read from kss
        #  msh and ELE: set from the parameter symmetric (readmsh)
        #   eltypes automatically converted to PMD according to 
        #   self.symmetric
        #  NOD: isn't set
        #writing:
        #  I1: written to kss
        #  msh: eltypes converted to MSH format according to symmetric
        #  NOD,ELE: isn't used
        self.crit = 3.0
        #reading:
        #  I1: read from crit 
        #  msh,NOD,ELE: isn't set 
        #writing:
        #  I1: written to crit
        #  msh,NOD,ELE: isn't used
    def readmsh(self,filename,b=None,symmetric=False,associateelements=True):
        """Reads mesh from filename (*.msh).

        it will read the physical entity "mshmodelnum" (default 100),
        which must contain every node (which will be used) and
        every element.

        Optional parameter b of type bound will be filled with
        all other physical (!=mshmodelnum) nodes.

        example:
            gmsh exports these physical entities: 100,1,2,101,200

            then readmsh will read 100 into self.elements and
            self.nodes, and optionally fills "b" with entities
            1,2,101 and 200. The assumption is, that entity 100
            will contain every node and element used in entities
            1,2,101 and 200. Also the nodes and elements in 100 
            must be consecutively sorted (use mshsort for this
            purpose)

        it will convert msh types to PMD types, so that
        self.elements only contains PMD types

        symmetric.... is the problem rot symmetric?
        if yes, readmsh() will automatically convert triangles
        to rottriangles and quadrangles to rotquadrangles.
        """
        self.clean()
        f=file(filename,"r")
        l=f.readline()
        check(l,"$NOD\n")
        l=f.readline()
        nnod=int(l)
        l=f.readline()
        n=1
        M=1
        xl=+M;xu=-M;yl=M;yu=-M;zl=M;zu=-M;
        self.symmetric=symmetric
        while l:
            x=string.split(l)
            p=[float(a) for a in x]
            if p[0] != n: 
                error("node-number mismatch (n=%d;p[0]=%d)"\
                    %(n,p[0]),2)
            if p[1]<xl:xl=p[1]
            if p[1]>xu:xu=p[1]
            if p[2]<yl:yl=p[2]
            if p[2]>yu:yu=p[2]
            if p[3]<zl:zl=p[3]
            if p[3]>zu:zu=p[3]
            if x[3]!="0": self.is2d=False
            self.nodes.append((int(p[0]),p[1],p[2],p[3]))
            l=f.readline()
            if n==nnod: break
            n+=1
        if symmetric and not self.is2d:
            error("symmetric and it isn't 2D!",2)
        if b!=None:
            b.is2d=self.is2d
        check(l,"$ENDNOD\n")
        l=f.readline()
        check(l,"$ELM\n")
        l=f.readline()
        nelm=int(l)
        l=f.readline()
        n=1
        pmdelm=0
        faces=[]
        while l and nelm != 0:
            x=string.split(l)
            p=[int(a) for a in x]
            if p[0] != n: error("elm-number mismatch",2)
            if p[2] == mshmodelnum:
                if p[1] == mshtriangle: 
                    if not self.is2d:
                        error("2D element in 3D mesh",2)
                    pmdelm+=1
                    if symmetric:
                        eltype=pmdtrianglerot
                    else:
                        eltype=pmdtriangle
                    self.elements.append((pmdelm,eltype,
                        p[5],p[6],p[7]))
                elif p[1] == mshtriangle2:
                    if not self.is2d:
                        error("2D element in 3D mesh",2)
                    pmdelm+=1
                    if symmetric:
                        eltype=pmdtrianglerot
                    else:
                        eltype=pmdtriangle
                    self.elements.append((pmdelm,eltype,
                        p[5],p[6],p[7],p[8],p[9],p[10]))
                elif p[1] == mshquadrangle:
                    if not self.is2d:
                        error("2D element in 3D mesh",2)
                    pmdelm+=1
                    if symmetric:
                        eltype=pmdquadranglerot
                    else:
                        eltype=pmdquadrangle
                    self.elements.append((pmdelm,eltype,
                        p[5],p[6],p[7],p[8]))
                elif p[1] == mshquadrangle2:
                    if not self.is2d:
                        error("2D element in 3D mesh",2)
                    pmdelm+=1
                    if symmetric:
                        eltype=pmdquadranglerot
                    else:
                        eltype=pmdquadrangle
                    self.elements.append((pmdelm,eltype,
                        p[5],p[6],p[7],p[8],
                        p[9],p[10],p[11],p[12]))
                    faces.append(p[13])
                elif p[1] == mshtetrahedron:
                    if self.is2d:
                        error("3D element in 2D mesh",2)
                    pmdelm+=1
                    eltype=pmdtetrahedron
                    self.elements.append((pmdelm,eltype,
                        p[5],p[6],p[7],p[8]))
                elif p[1] == mshtetrahedron2:
                    if self.is2d:
                        error("3D element in 2D mesh",2)
                    pmdelm+=1
                    eltype=pmdtetrahedron
                    self.elements.append((pmdelm,eltype,
                        p[5],p[6],p[7],p[8],
                        p[9],p[10],p[11],p[12],p[13],p[14]))
                elif p[1] == mshhexahedron:
                    if self.is2d:
                        error("3D element in 2D mesh",2)
                    pmdelm+=1
                    eltype=pmdhexahedron
                    el=[pmdelm,eltype]
                    el.extend(p[5:5+8])
                    self.elements.append(tuple(el))
                elif p[1] == mshprism:
                    if self.is2d:
                        error("3D element in 2D mesh",2)
                    pmdelm+=1
                    eltype=pmdprism
                    el=[pmdelm,eltype]
                    el.extend(p[5:5+6])
                    self.elements.append(tuple(el))
                else:
                    error("unsupported el %d"%(p[1]),3)
            elif p[2] < mshmodelnum:
                if b!=None:
                    b.handle2(p)
            else:
                if b!=None:
                    b.handleelement(p)
            l=f.readline()
            #if n==nelm: break
            if l[0]=="$": break
            n+=1
        check(l,"$ENDELM\n")
        l=f.readline()
        if l != "": error("extra lines at the end of file",2)
        f.close()
        eps=0.001
        self.boundbox=(xl-eps, xu+eps, yl-eps, yu+eps, zl-eps, zu+eps)
        self.removecentralnodes(faces)
        if b!=None and associateelements:
            b.associateelements(self.elements)
    def readmsh2(self,filename,b=None,symmetric=False,associateelements=True):
        """Reads mesh from filename (*.msh). Version 2.0

        it will read the physical entity "mshmodelnum" (default 100),
        which must contain every node (which will be used) and
        every element.

        Optional parameter b of type bound will be filled with
        all other physical (!=mshmodelnum) nodes.

        example:
            gmsh exports these physical entities: 100,1,2,101,200

            then readmsh will read 100 into self.elements and
            self.nodes, and optionally fills "b" with entities
            1,2,101 and 200. The assumption is, that entity 100
            will contain every node and element used in entities
            1,2,101 and 200. Also the nodes and elements in 100 
            must be consecutively sorted (use mshsort for this
            purpose)

        it will convert msh types to PMD types, so that
        self.elements only contains PMD types

        symmetric.... is the problem rot symmetric?
        if yes, readmsh() will automatically convert triangles
        to rottriangles and quadrangles to rotquadrangles.
        """
        self.clean()
        f=file(filename,"r")
        l=f.readline()
        l=f.readline()
        l=f.readline()
        l=f.readline()
        check(l,"$Nodes\n")
        l=f.readline()
        nnod=int(l)
        l=f.readline()
        n=1
        M=1
        xl=+M;xu=-M;yl=M;yu=-M;zl=M;zu=-M;
        self.symmetric=symmetric
        while l:
            x=string.split(l)
            p=[float(a) for a in x]
            if p[0] != n: 
                error("node-number mismatch (n=%d;p[0]=%d)"\
                    %(n,p[0]),2)
            if p[1]<xl:xl=p[1]
            if p[1]>xu:xu=p[1]
            if p[2]<yl:yl=p[2]
            if p[2]>yu:yu=p[2]
            if p[3]<zl:zl=p[3]
            if p[3]>zu:zu=p[3]
            if x[3]!="0": self.is2d=False
            self.nodes.append((int(p[0]),p[1],p[2],p[3]))
            l=f.readline()
            if n==nnod: break
            n+=1
        if symmetric and not self.is2d:
            error("symmetric and it isn't 2D!",2)
        if b!=None:
            b.is2d=self.is2d
        check(l,"$EndNodes\n")
        l=f.readline()
        check(l,"$Elements\n")
        l=f.readline()
        nelm=int(l)
        l=f.readline()
        n=1
        pmdelm=0
        faces=[]
        while l and nelm != 0:
            x=string.split(l)
            p=[int(a) for a in x]
            if p[0] != n: error("elm-number mismatch",2)
            if p[2] == mshmodelnum:
                if p[1] == mshtriangle: 
                    if not self.is2d:
                        error("2D element in 3D mesh",2)
                    pmdelm+=1
                    if symmetric:
                        eltype=pmdtrianglerot
                    else:
                        eltype=pmdtriangle
                    self.elements.append((pmdelm,eltype,
                        p[5],p[6],p[7]))
                elif p[1] == mshtriangle2:
                    if not self.is2d:
                        error("2D element in 3D mesh",2)
                    pmdelm+=1
                    if symmetric:
                        eltype=pmdtrianglerot
                    else:
                        eltype=pmdtriangle
                    self.elements.append((pmdelm,eltype,
                        p[5],p[6],p[7],p[8],p[9],p[10]))
                elif p[1] == mshquadrangle:
                    if not self.is2d:
                        error("2D element in 3D mesh",2)
                    pmdelm+=1
                    if symmetric:
                        eltype=pmdquadranglerot
                    else:
                        eltype=pmdquadrangle
                    self.elements.append((pmdelm,eltype)+tuple(p[3:]))
                elif p[1] == mshquadrangle2:
                    if not self.is2d:
                        error("2D element in 3D mesh",2)
                    pmdelm+=1
                    if symmetric:
                        eltype=pmdquadranglerot
                    else:
                        eltype=pmdquadrangle
                    self.elements.append((pmdelm,eltype,
                        p[5],p[6],p[7],p[8],
                        p[9],p[10],p[11],p[12]))
                    faces.append(p[13])
                elif p[1] == mshtetrahedron:
                    if self.is2d:
                        error("3D element in 2D mesh",2)
                    pmdelm+=1
                    eltype=pmdtetrahedron
                    self.elements.append((pmdelm,eltype,
                        p[5],p[6],p[7],p[8]))
                elif p[1] == mshhexahedron:
                    if self.is2d:
                        error("3D element in 2D mesh",2)
                    pmdelm+=1
                    eltype=pmdhexahedron
                    el=[pmdelm,eltype]
                    el.extend(p[5:5+8])
                    self.elements.append(tuple(el))
                elif p[1] == mshprism:
                    if self.is2d:
                        error("3D element in 2D mesh",2)
                    pmdelm+=1
                    eltype=pmdprism
                    el=[pmdelm,eltype]
                    el.extend(p[5:5+6])
                    self.elements.append(tuple(el))
                else:
                    error("unsupported el %d"%(p[1]),3)
            elif p[2] < mshmodelnum:
                if b!=None:
                    b.handle2(p)
            else:
                if b!=None:
                    b.handleelement(p)
            l=f.readline()
            if n==nelm: break
            if l[0]=="$": break
            n+=1
        l=f.readline()
        check(l,"$EndElements\n")
        l=f.readline()
        if l != "": error("extra lines at the end of file",2)
        f.close()
        eps=0.001
        self.boundbox=(xl-eps, xu+eps, yl-eps, yu+eps, zl-eps, zu+eps)
        self.removecentralnodes(faces)
        if b!=None and associateelements:
            b.associateelements(self.elements)
    def writexda(self,filename,verbose=True,b=None):
        """Writes mesh to filename (*.xda).

        We try to be byte to byte compatible with the xda output from libmesh
        (so I use the same tabs and spaces as libmesh does).

        """
        up=progressbar.MyBar("Writing mesh to %s:"%filename)
        if verbose: up.init(len(self.nodes)+2*len(self.elements))
        mapping=[]
        blocks={}
        sew=0
        c=0
        for e in self.elements:
            p=[n-1 for n in e[2:]]
            t=e[1]
            if t==pmdtriangle:
                elt=libmeshtriangle
            elif t==pmdquadrangle:
                elt=libmeshquadrangle
            elif t==pmdtetrahedron:
                if len(p)==4:
                    elt=libmeshtetrahedron
                else:
                    assert len(p)==10
                    elt=libmeshtetrahedron_quad
            elif t==pmdhexahedron:
                elt=libmeshhexahedron
            elif t==pmdprism:
                elt=libmeshprism
            else:
                error("Unimplemented yet t=%d."%(t),2)
            if elt in blocks:
                blocks[elt].append(p)
            else:
                blocks[elt]=[p]
            mapping.append((elt,len(blocks[elt])))
            sew+=len(p)
            c+=1
            if verbose: up.update(c)
        nels=[len(blocks[t]) for t in blocks.keys()]
        map2=[]
        for x in mapping:
            k=[y for y in blocks.keys() if y<x[0]]
            n=0
            for t in k:
                n+=len(blocks[t])
            map2.append(n+x[1]-1)
        mapping=map2
        #print mapping

        bs=[]
        if b:
            for key in b.f.keys(): 
                if key >=mshmodelnum: continue
                bo=[ [mapping[el-1],side-1,key] for (el,side) in 
                    b.findelements(key,self.elements)]
                bs.extend(bo)

        f=file(filename,"w")
        f.write("DEAL 003:003\n")
        f.write("%d  # Num. Elements\n"%len(self.elements))
        f.write("%d  # Num. Nodes\n"%len(self.nodes))
        f.write("%d  # Sum of Element Weights\n"%(sew))
        f.write("%d  # Num. Boundary Conds.\n"%(len(bs)))
        f.write("%d  # String Size (ignore)\n"%(65536))
        f.write("%d  # Num. Element Blocks.\n"%len(blocks))
        f.write(("%d "*len(blocks.keys()))%tuple(blocks.keys())+
            "    # Element types in each block.\n")
        f.write(("%d "*len(nels))%tuple(nels)+
            "    # Num. of elements in each block.\n")
        f.write("Id String\n")
        f.write("Title String\n")
        for block in blocks.values():
            for el in block:
                f.write(("%d "*len(el))%tuple(el)+"\n")
                c+=1
                if verbose: up.update(c)
        for node in self.nodes:
                f.write("%e     %e  %e  \n"%tuple(node[1:]))
                c+=1
                if verbose: up.update(c)
        for line in bs:
            f.write(("%d "*len(line))%tuple(line)+"\n")
    def writemsh(self,filename,verbose=True):
        """Writes mesh to filename (*.msh).

        """
        up=progressbar.MyBar("Writing mesh to %s:"%filename)
        if verbose: up.init(len(self.nodes)+len(self.elements))
        f=file(filename,"w")
        l=f.write("$NOD\n")
        l=f.write("%d\n"%len(self.nodes))
        for node in self.nodes:
            if self.is2d:
                f.write("%d %f %f %d\n"%node)
            else:
                f.write("%d %f %f %f\n"%node)
            if verbose: up.update(node[0])
        l=f.write("$ENDNOD\n")
        l=f.write("$ELM\n")
        l=f.write("%d\n"%len(self.elements))
        for el in self.elements:
            if el[1]==pmdtriangle or el[1]==pmdtrianglerot:
                eltype=mshtriangle
                number_of_nodes=3
                #if we want 2nd order, fix it here
            elif el[1]==pmdquadrangle or el[1]==pmdquadranglerot:
                eltype=mshquadrangle
                number_of_nodes=4
                #if we want 2nd order, fix it here
            elif el[1]==pmdsemiloof:
                #eltype=mshtriangle
                eltype=mshquadrangle
                number_of_nodes=4
            elif el[1]==pmdhexahedron:
                eltype=mshhexahedron
                number_of_nodes=8
                #if we want 2nd order, fix it here
            elif el[1]==pmdtetrahedron:
                if len(el[2:])==4:
                    eltype=mshtetrahedron
                    number_of_nodes=4
                elif len(el[2:])==10:
                    eltype=mshtetrahedron2
                    number_of_nodes=10
                else:
                    assert False
            elif el[1]==pmdprism:
                eltype=mshprism
                number_of_nodes=6
                #if we want 2nd order, fix it here
            else:
                error("unsupported eltype type=%s"\
                    %(repr(el[1])),3)
            n=[ #elm-number 
                el[0],
                #elm-type 
                eltype,
                #reg-phys
                mshmodelnum,
                #reg-elem 
                0,
                #number-of-nodes
                number_of_nodes]
                #node-number-list
            n.extend(el[2:2+number_of_nodes])
            f.write("%d "*len(n)%tuple(n))
            f.write("\n")
            if verbose: up.update(len(self.nodes)+el[0])
        l=f.write("$ENDELM\n")
        f.close()
    def writemsh2(self,filename):
        """Writes mesh to filename (*.msh). Version 2.0

        """
        f=file(filename,"w")
        l=f.write("$MeshFormat\n")
        l=f.write("2.0 0 8\n")
        l=f.write("$EndMeshFormat\n")
        l=f.write("$Nodes\n")
        l=f.write("%d\n"%len(self.nodes))
        for node in self.nodes:
            if self.is2d:
                f.write("%d %f %f %d\n"%node)
            else:
                f.write("%d %f %f %f\n"%node)
        l=f.write("$EndNodes\n")
        l=f.write("$Elements\n")
        l=f.write("%d\n"%len(self.elements))
        for el in self.elements:
            if el[1]==pmdtriangle or el[1]==pmdtrianglerot:
                eltype=mshtriangle
                number_of_nodes=3
                #if we want 2nd order, fix it here
            elif el[1]==pmdquadrangle or el[1]==pmdquadranglerot:
                eltype=mshquadrangle
                number_of_nodes=4
                #if we want 2nd order, fix it here
            elif el[1]==pmdsemiloof:
                #eltype=mshtriangle
                eltype=mshquadrangle
                number_of_nodes=4
            elif el[1]==pmdhexahedron:
                eltype=mshhexahedron
                number_of_nodes=8
                #if we want 2nd order, fix it here
            elif el[1]==pmdtetrahedron:
                eltype=mshtetrahedron
                number_of_nodes=4
                #if we want 2nd order, fix it here
            elif el[1]==pmdprism:
                eltype=mshprism
                number_of_nodes=6
                #if we want 2nd order, fix it here
            else:
                error("unsupported eltype type=%s"\
                    %(repr(el[1])),3)
            n=[ #elm-number 
                el[0],
                #elm-type 
                eltype,
                #reg-phys
                mshmodelnum,
                #reg-elem 
                0,
                #number-of-nodes
                number_of_nodes]
                #node-number-list
            n.extend(el[2:2+number_of_nodes])
            f.write("%d "*len(n)%tuple(n))
            f.write("\n")
        l=f.write("$EndElements\n")
        f.close()
    def readNOD(self,filename,scale=1.0):
        """Read nodes from filename (*.NOD). 
        """
        f=file(filename)
        l=f.readline()
        self.nodes=[]
        self.is2d=True
        while l:
            x=string.split(l)
            node=(int(x[0]),
                myfloat(x[1])*scale,
                myfloat(x[2])*scale,
                myfloat(x[3])*scale)
            if node[3]!=0.0:
                self.is2d=False
            self.nodes.append(node)
            l=f.readline()
    def writeNOD(self,filename):
        """Write nodes to filename (*.NOD). 
        """
        f=file(filename,"w")
        for nod in self.nodes:
            f.write("%d %f %f %f\n"%tuple(nod))
    def readELE(self,filename,symmetric=False):
        """Read elements from filename (*.ELE). 
        """
        #format:first line don't know yet 
        #(n,NDIM,nnod,n1,n2,...,n_nnod,T1,T2,...,T_nnod)...
        #n....... element number
        #NDIM ... probably dimension of the problem 2 or 3
        #nnod ... number of nodes of the element
        #n1...n_nnod ... nodes
        #T1...T_nnod ... THICK (probably only for 2D problems)
        f=file(filename)
        l=f.readline()
        data=[]
        for l in f.readlines(): data.extend(string.split(l))
        self.elements=[]
        n=1
        pos=0
        while pos<len(data):
            if n!=int(data[pos]):
                error("element number mischmatch",2)
            nnod=int(data[pos+2])
            x=[int(i) for i in data[pos:pos+3+nnod]]
            ndim=x[1]
            if nnod==6:
                if symmetric:
                    ite=pmdtrianglerot
                else:
                    ite=pmdtriangle
            elif nnod==8:
                if symmetric:
                    ite=pmdquadranglerot
                else:
                    ite=pmdquadrangle
            else:
                error("unsupported element",2)
            el=[n,ite]+x[3:3+nnod]
            self.elements.append(el)
            n+=1
            if ndim==2:
                nnod*=2
            pos+=3+nnod

    def renumber_elements(self):
        self.writeELE("t1.ELE")
        if os.access("XELM.ELE",os.F_OK): os.remove("XELM.ELE")
        os.spawnv(os.P_WAIT,"/home/ondra/pmd/PMD/xelm",["xelm","t1.ELE"])
        print "readELE"
        self.readELE("XELM.ELE")
        print "removing"
        os.remove("XELM.ELE")
        os.remove("t1.ELE")
        el=[]
        for e in self.elements:
            el.append(e[:-3])
        self.elements=el

    def readELE2(self,filename):
        """Read elements from filename (*.ELE). 
        """
        #format:first line don't know yet 
        #(n,nnod-1,n1,n2,...,n_nnod)...
        #n....... element number
        #nnod ... number of nodes of the element
        #n1...n_nnod ... nodes
        f=file(filename)
        data=[]
        for l in f.readlines(): data.extend(string.split(l))
        self.elements=[]
        n=1
        while data:
        #   if n!=int(data[0]):
        #       error("element number mischmatch",2)
            nnod=int(data[1])+1
            x=[int(i) for i in data[0:2+nnod]]
            if nnod==6:
                    ite=pmdtriangle
            elif nnod==8:
                    ite=pmdhexahedron
            else:
                error("unsupported element",2)
            el=[n,ite]+x[2:2+nnod]
            self.elements.append(el)
            n+=1
            data=data[2+nnod:]
    def sortnodes(self):
        def mycmp(a,b):
            if a[0] > b[0]:
                return 1
            elif a[0] < b[0]:
                return -1
            else:
                return 0
        self.nodes.sort(mycmp)
    def writeELE(self,filename):
        """Write nodes to filename (*.ELE). 
        """
        f=file(filename,"w")
        f.write("%d %d %d %d\n"%(0,2,8,2))
        if self.is2d:
            ndim=2
        else:
            ndim=3
        for el in self.elements:
            nods=tuple(el[2:])+(0,0,0)
            nums=(el[0],ndim,len(nods))+tuple(nods)
            f.write((" %d"*len(nums)+"\n")%nums)
            if self.is2d:
                nums=[1.0]*len(nods)
                f.write((" %.2f"*len(nums)+"\n")%tuple(nums))
    def readpmd(self,filename):
        """Read the mesh from filename (*.I1). 

        Can read I1, which was produced by writepmd(). 

        Should be able to read some hand written I1s from PMD Example
        Manual. Sometimes you will have to tweak the parser to read
        your syntax. 

        I suggest to only use the writepmd() I1 syntax. This
        is easily readable/writable. 
        
        """
        self.clean()
        f=file(filename,"r")
        def C(x):
            if x[0]==";":
                x=f.readline()
            return x
        x=string.split(C(f.readline()))
        check(x[0],"IP")
        p=[int(a) for a in x[1:]]
        nelements=p[0]; nnodes=p[1]; ited=p[2]; kss=p[9]
        if kss==2:
            self.symmetric=True
        else:
            self.symmetric=False

        x=string.split(C(f.readline()))
        check(x[0],"RP")
        p=[float(a) for a in x[1:]]
        self.crit=p[0]; scale=p[1]; thdef=p[2]; self.boundbox=p[3:9]

        x=string.split(C(f.readline()))
        check(x[0],"XY")
        if len(x)>1: #old format....
            dict={}
            def ev1(_str,loc,toks):
                out= range(int(toks[0]),int(toks[2])+1)
                out=[str(a) for a in out]
            #   print toks, "->", out
                return out
            def ev2(str,loc,toks):
                out=int(toks[0])*toks.asList()[2:]
            #   print toks, "->", out
                return out
            def ev3(str,loc,toks):
                key=toks[1] 
                out=toks.asList()[2:len(toks)-2]
                d=[]
                for a in out:
                    num=eval("%s"%(a))
                    d.append(num)
                dict[key]=d
                #print toks, "->", out
                return out
            def ev4(_str,loc,toks):
                if toks[1][0] in nums:
                    n=eval("%s"%(toks[1]))
                    key=toks[2] 
                else:
                    n=0
                    key=toks[1] 
                assert(dict.has_key(key))
                #increment dict[key] by n
                out=[str(a+n) for a in dict[key]]
                #print toks, "dict->", out
                return out
            point = Literal(".")
            e     = CaselessLiteral("E")
            inum  = Combine(Word("+-"+nums,nums)+ 
                    Optional(e+Word( "+-"+nums, nums ) ) )
            fnum  = Combine( Word( "+-"+nums, nums ) + 
                    Optional(point+Optional(Word(nums))) +
                    Optional(e+Word("+-"+nums,nums ) ) )
            semi  = Literal(";")
            lpar  = Literal("(").suppress()
            rpar  = Literal(")").suppress()
            callpar = Literal("=")+Optional(inum)+Word(alphas,max=1)
            defpar  = Literal("=")+Word(alphas,max=1)
            comment = semi + Optional(restOfLine)

            terms=Forward()
            atom=fnum | (lpar+ terms+rpar)
            seq=(atom+":"+atom).setParseAction(ev1)
            rep=(inum+"*"+atom).setParseAction(ev2)
            terms << OneOrMore(rep | seq | atom)
            numlist=OneOrMore(terms|(defpar+terms+
                defpar).setParseAction(ev3) | 
                callpar.setParseAction(ev4))

            nodX = Group("X"+numlist)
            #nodX = Group("X"+numlist+LineEnd())
            nodY = Group("Y"+numlist)
            nodZ = Group("Z"+numlist)
            nodes=Group(Dict(Literal("XY")+Group("N"+numlist)+nodX+nodY+Optional(nodZ)))
            element= Group(Literal("EL")+Optional("T"+inum)+"E"+numlist+"N"+
                OneOrMore(numlist))
            CN=Group("CN"+restOfLine)
            elements=Group(OneOrMore(element)).setResultsName("ELs")

            EN=Literal("EN")
            end=Group(EN+EN)
            grammar=Dict(nodes+elements+end+StringEnd())
            grammar.ignore(comment)
            grammar.ignore(CN)

            data=[string.join(x)+"\n"]
            data.extend(f.readlines())
            data=string.join(data)
            tokens=""
            try:
                tokens=grammar.parseString(data)
            except ParseException, err:
                error("\n"+err.line+"\n"+" "*(err.column-1)+\
                    "^\n" + repr(err),2)

            if "Z" in tokens["XY"].keys():
                Zn=map(float,tokens["XY"]["Z"])
                self.is2d=False
            else:
                Zn=[0.0]*len(tokens["XY"]["X"])

            self.nodes=zip(range(1,len(Zn)+1),
                map(float,tokens["XY"]["X"]),
                map(float,tokens["XY"]["Y"]),
                Zn)

            elm={}
            els= tokens["ELs"]
            for el in els:
                i=1
                if el[i]=="T":
                    type=int(el[i+1])
                    i+=2
                else:
                    type=ited
                if type==pmdhexahedron:
                    nnum=20
                    nskip=20
                elif type==pmdtriangle:
                    nnum=3
                    nskip=6
                elif type==pmdquadrangle:
                    nnum=4
                    nskip=8
                elif type==pmdsemiloof:
                    nnum=8
                    nskip=8
                else:
                    error("unsupported type. "\
                        "type=%d"%(type),3)
                el=el.asList()[i:]
                n=1
                while el[n]!="N": n+=1
                n+=1
                for e in el[1:n-1]:
                    elm[int(e)]=(int(e),type)+tuple(
                        map(int,el[n:n+nnum]))
                    n+=nskip
            n=1
            for i,e in elm.iteritems():
                self.elements.append(e)
                assert(i==n)
                n+=1
        else:
            for n in range(1,nnodes+1):
                x=string.split(C(f.readline()))
                check(x[0],"C")
                if int(x[1])!=n: error("node number mismatch",2)
                if len(x)==4:
                    node=(n,float(x[2]),float(x[3]),
                        float(0))
                else:
                    node=(n,float(x[2]),float(x[3]),
                        float(x[4]))
                    self.is2d=False
                self.nodes.append(node)

            x=string.split(C(f.readline()))
            for n in range(1,nelements+1):
            #   print x
                check(x[0],"EL")
                check(x[1],"T")
                T=int(x[2])
                check(x[3],"E")
                if int(x[4])!=n: error("node number mismatch",2)
                check(x[5],"N")
                el=[n,T]
                el.extend([int(a) for a in x[6:]])
                x=string.split(C(f.readline()))
                if x[0] != "EL" and x[0] != "EN":
                    el.extend([int(a) for a in x])
                    x=string.split(C(f.readline()))
                self.elements.append(el)
            check(x[0],"EN")
            x=string.split(C(f.readline()))
            check(x[0],"EN")
            f.close()

        #apply scale factor
        self.nodes=[(i[0],i[1]*scale,i[2]*scale,i[3]*scale) for i in
            self.nodes]
        self.boundbox=[i*scale for i in self.boundbox]
    def writepmd(self,filename):
        "Writes the mesh to filename (*.I1)."
        f=file(filename,"w")
        ited=4;
        if self.symmetric:
            kss=2
        else:
            kss=-1
        f.write("IP %d %d %d 0 0 0 0 0 0 %d\n"
            %(len(self.elements),len(self.nodes),ited,kss));
        scale=1.0;thdef=1;
        str=[self.crit,scale,thdef]
        str.extend(self.boundbox)
        f.write("RP %.2f %.4f %.3f %.3f %.3f %.3f %.3f %.3f %.3f 0\n"
            %tuple(str))
        f.write("XY\n");
        for p in self.nodes:
            if self.is2d:
                f.write("  C %d %.17f %.17f\n"%(p[0],p[1],p[2]))
            else:
                f.write("  C %d %.17f %.17f %.17f\n"%
                    (p[0],p[1],p[2],p[3]))
        for p in self.elements:
            f.write("EL T %d E %d N"%(p[1],p[0]))
            for a in p[2:]: f.write(" %d"%(a))
            f.write("\n")
        f.write("EN\n");
        f.write("EN\n");
        f.close()
    def readxt2sSTR(self,filename):
        "Reads temperature data from filename (*.STR) to scalars."
        f=file(filename,"r")
        l=f.readline()
        l=f.readline()
        check(l,"$TEMPERATURE\n")
        l=f.readline()
        n=0
        self.scalars=[]
        while l:
            n+=1
            x=string.split(l)
            p=[float(a) for a in x]
            if p[0] != n: error("node-number mismatch",2)
            self.scalars.append(p[1])
            if p[2] != 0: error("2nd number is not zero",2)
            if p[3] != 0: error("2nd number is not zero",2)
            l=f.readline()
        f.close()
    def readstr2STR(self,filename):
        "Reads temperature data from filename (*.STR) to scalars."
        f=file(filename,"r")
        l=f.readline()
        l=f.readline()
        check(l,"$DISPLACEMENT\n")
        l=f.readline()
        n=0
        self.vectors=[]
        while l:
            n+=1
            x=string.split(l)
            p=[float(a) for a in x]
            if p[0] != n: error("node-number mismatch",2)
            self.vectors.append(p[1:4])
            l=f.readline()
            if n==len(self.nodes): break
        check(l,"$ELEMENT\n")
        l=f.readline()
        n=0
        self.elementdata=[]
        while l:
            n+=1
            x=string.split(l)
            check(x[0],"*IE")
            if int(x[1]) != n: error("node-number mismatch",2)
            check(int(x[2]),1)
            check(float(x[3]),0)
            l=f.readline()
            check(l,"*STRESS\n")
            for i in range(24):
                l=f.readline()
            l=f.readline()
            check(l,"*HMH\n")
            nums=[]
            for i in range(4):
                l=f.readline()
                x=string.split(l)
                p=[float(a) for a in x]
                nums.extend(p)
            self.elementdata.append(nums)
            l=f.readline()
            check(l,"*TAU\n")
            for i in range(4):
                l=f.readline()
            l=f.readline()
            check(l,"*SCALAR1\n")
            for i in range(4):
                l=f.readline()
            l=f.readline()
            check(l,"*SCALAR2\n")
            for i in range(4):
                l=f.readline()
            l=f.readline()
        f.close()
    def readstr3STR(self,filename):
        "Reads temperature data from filename (*.STR) to scalars."
        f=file(filename,"r")
        l=f.readline()
        l=f.readline()
        check(l,"$DISPLACEMENT\n")
        l=f.readline()
        n=0
        self.vectors=[]
        while l:
            n+=1
            x=string.split(l)
            p=[float(a) for a in x]
            if p[0] != n: error("node-number mismatch",2)
            self.vectors.append(p[1:4])
            l=f.readline()
            if n==len(self.nodes): break
        check(l,"$ELEMENT\n")
        l=f.readline()
        n=0
        self.elementdata=[]
        while l:
            n+=1
            x=string.split(l)
            check(x[0],"*IE")
            if int(x[1]) != n: error("node-number mismatch",2)
            check(int(x[2]),1)
            check(float(x[3]),0)
            l=f.readline()
            check(l,"*STRESS\n")
            for i in range(24):
                l=f.readline()
            l=f.readline()
            check(l,"*HMH\n")
            nums=[]
            for i in range(4):
                l=f.readline()
                x=string.split(l)
                p=[float(a) for a in x]
                nums.extend(p)
            self.elementdata.append(nums)
            l=f.readline()
            check(l,"*TAU\n")
            for i in range(4):
                l=f.readline()
            l=f.readline()
        f.close()
    def writevectorspos(self,filename,infotext="PMD_vectors"):
        """Writes scalars together with nodes and elements to filename.

        Optional parameter infotext specifies the name of the view 
        in the pos file.

        1) Associates a scalar with every node, so gmsh
        shows points.
        2) Associates a scalar with every node of all elements, so gmsh
        fills the whole element (triangle, quadrangle etc.) with an
        extrapolated color.

        You can set visibility in gmsh (only points, only triangles...).
        """
        if len(self.nodes) != len(self.vectors):
            error("Different number of nodes and vectors!")
        f=file(filename,"w")
        f.write("$PostFormat\n")
        #1.3 file-type data-size
        f.write("%g %d %d\n"%(1.2,0,8))
        f.write("$EndPostFormat\n")
        f.write("$View\n")
        f.write("%s %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d %d\n"%(
    #view-name nb-time-steps
            infotext, 1, 
    #nb-scalar-points nb-vector-points nb-tensor-points
            0,len(self.nodes),0, 
    #nb-scalar-lines nb-vector-lines nb-tensor-lines
            0,0,0,
    #nb-scalar-triangles nb-vector-triangles nb-tensor-triangles
            0,0,0,
            #0,self.getnumtriangles(),0,
    #nb-scalar-quadrangles nb-vector-quadrangles nb-tensor-quadrangles
            0,0,0,
            #0,self.getnumquadrangles(),0,
    #nb-scalar-tetrahedra nb-vector-tetrahedra nb-tensor-tetrahedra
            0,0,0,
    #nb-scalar-hexahedra nb-vector-hexahedra nb-tensor-hexahedra
            0,0,0,
    #nb-scalar-prisms nb-vector-prisms nb-tensor-prisms
            0,0,0,
    #nb-scalar-pyramids nb-vector-pyramids nb-tensor-pyramids
            0,0,0,
    #nb-text2d nb-text2d-chars nb-text3d nb-text3d-chars
            0,0,0,0
            ))
        #time-step-values
        f.write("%d\n"%(0))
        #< scalar-point-value > ...
        for node in self.nodes:
            n=(self.getxyz(node[0]),)
            T=(self.getvector(node[0]),)
            f.write(formatpos2(n,T))
        f.write("$EndView\n")
        f.close()
    def writevectorspos3(self,filename,vectorfield,infotext="PMD_vectors"):
        self.vectors=vectorfield
        self.writevectorspos(filename,infotext)
    def writescalarspos3(self,filename,scalars,infotext="PMD_scalars"):
        self.scalars=scalars
        self.writescalarspos(filename,infotext)
    def writescalarspos(self,filename,infotext="PMD_scalars"):
        """Writes self.scalars to *.pos.

        Optional parameter infotext specifies the name of the view 
        in the pos file.

        1) Associates a scalar with every node, so gmsh shows points.
        2) Associates a scalar with every node of all elements, so gmsh fills
           the whole element (triangle, quadrangle, tetrahedra, ...  etc.)
           with an extrapolated color.

        You can set visibility in gmsh (only points, only triangles...).
        """
        if len(self.nodes) != len(self.scalars):
            error("Different number of nodes and scalars!")
        up=progressbar.MyBar("Writing scalar field to %s:"%filename)
        up.init(len(self.nodes)+len(self.elements))
        f=file(filename,"w")
        f.write("$PostFormat\n")
        #1.3 file-type data-size
        f.write("%g %d %d\n"%(1.2,0,8))
        f.write("$EndPostFormat\n")
        f.write("$View\n")
        f.write("%s %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d %d\n"%(
    #view-name nb-time-steps
            infotext, 1, 
    #nb-scalar-points nb-vector-points nb-tensor-points
            len(self.nodes),0,0, 
    #nb-scalar-lines nb-vector-lines nb-tensor-lines
            0,0,0,
    #nb-scalar-triangles nb-vector-triangles nb-tensor-triangles
            self.getnumtriangles(),0,0,
    #nb-scalar-quadrangles nb-vector-quadrangles nb-tensor-quadrangles
            self.getnumquadrangles(),0,0,
    #nb-scalar-tetrahedra nb-vector-tetrahedra nb-tensor-tetrahedra
            self.getnumtetrahedra(),0,0,
    #nb-scalar-hexahedra nb-vector-hexahedra nb-tensor-hexahedra
            self.getnumhexahedra(),0,0,
    #nb-scalar-prisms nb-vector-prisms nb-tensor-prisms
            self.getnumprisms(),0,0,
    #nb-scalar-pyramids nb-vector-pyramids nb-tensor-pyramids
            0,0,0,
    #nb-text2d nb-text2d-chars nb-text3d nb-text3d-chars
            0,0,0,0
            ))
        #time-step-values
        f.write("%d\n"%(0))
        c=0
        #< scalar-point-value > ...
        for node in self.nodes:
            n=(self.getxyz(node[0]),)
            T=(self.getscalar(node[0]),)
            f.write(formatpos(n,T))
            c+=1
            up.update(c)
        #< scalar-triangles-value > ...
        for el in self.elements:
            if el[1] in [pmdtriangle,pmdtrianglerot]: 
                n=(self.getxyz(el[2]),
                self.getxyz(el[3]),
                self.getxyz(el[4]))
                T=(self.getscalar(el[2]),
                self.getscalar(el[3]),
                self.getscalar(el[4]))
                f.write(formatpos(n,T))
        #< scalar-quadrangles-value > ...
        for el in self.elements:
            if el[1] == pmdquadrangle: 
                n=(self.getxyz(el[2]),
                self.getxyz(el[3]),
                self.getxyz(el[4]),
                self.getxyz(el[5]))
                T=(self.getscalar(el[2]),
                self.getscalar(el[3]),
                self.getscalar(el[4]),
                self.getscalar(el[5]))
                f.write(formatpos(n,T))
        #< scalar-tetrahedra-value > ...
        for el in self.elements:
            if el[1] == pmdtetrahedron: 
                n=(self.getxyz(el[2]),
                self.getxyz(el[3]),
                self.getxyz(el[4]),
                self.getxyz(el[5]))
                T=(self.getscalar(el[2]),
                self.getscalar(el[3]),
                self.getscalar(el[4]),
                self.getscalar(el[5]))
                f.write(formatpos(n,T))
                c+=1
                up.update(c)
        #< scalar-hexahedra-value > ...
        for el in self.elements:
            if el[1] == pmdhexahedron: 
                #print el
                n=(self.getxyz(el[2]),
                self.getxyz(el[3]),
                self.getxyz(el[4]),
                self.getxyz(el[5]),
                self.getxyz(el[6]),
                self.getxyz(el[7]),
                self.getxyz(el[8]),
                self.getxyz(el[9]))
                T=(self.getscalar(el[2]),
                self.getscalar(el[3]),
                self.getscalar(el[4]),
                self.getscalar(el[5]),
                self.getscalar(el[6]),
                self.getscalar(el[7]),
                self.getscalar(el[8]),
                self.getscalar(el[9]))
                f.write(formatpos(n,T))
        #< scalar-prisms-value > ...
        for el in self.elements:
            if el[1] == pmdprism: 
                #print el
                n=(self.getxyz(el[2]),
                self.getxyz(el[3]),
                self.getxyz(el[4]),
                self.getxyz(el[5]),
                self.getxyz(el[6]),
                self.getxyz(el[7]))
                T=(self.getscalar(el[2]),
                self.getscalar(el[3]),
                self.getscalar(el[4]),
                self.getscalar(el[5]),
                self.getscalar(el[6]),
                self.getscalar(el[7]))
                f.write(formatpos(n,T))
        f.write("$EndView\n")
        f.close()
    def writescalarspos2(self,filename,scaltime,infotext="PMD_scalars",dt=1.0):
        """Writes self.scalars to *.pos.

        Optional parameter infotext specifies the name of the view 
        in the pos file.

        1) Associates a scalar with every node, so gmsh shows points.
        2) Associates a scalar with every node of all elements, so gmsh fills
           the whole element (triangle, quadrangle, tetrahedra, ...  etc.)
           with an extrapolated color.

        You can set visibility in gmsh (only points, only triangles...).
        """
        f=file(filename,"w")
        f.write("$PostFormat\n")
        #1.3 file-type data-size
        f.write("%g %d %d\n"%(1.2,0,8))
        f.write("$EndPostFormat\n")
        f.write("$View\n")
        f.write("%s %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d %d\n"%(
    #view-name nb-time-steps
            infotext, len(scaltime), 
    #nb-scalar-points nb-vector-points nb-tensor-points
            #len(self.nodes),0,0, 
            0,0,0, 
    #nb-scalar-lines nb-vector-lines nb-tensor-lines
            0,0,0,
    #nb-scalar-triangles nb-vector-triangles nb-tensor-triangles
            self.getnumtriangles(),0,0,
    #nb-scalar-quadrangles nb-vector-quadrangles nb-tensor-quadrangles
            #self.getnumquadrangles(),0,0,
            0,0,0,
    #nb-scalar-tetrahedra nb-vector-tetrahedra nb-tensor-tetrahedra
            #self.getnumtetrahedra(),0,0,
            0,0,0,
    #nb-scalar-hexahedra nb-vector-hexahedra nb-tensor-hexahedra
            #self.getnumhexahedra(),0,0,
            0,0,0,
    #nb-scalar-prisms nb-vector-prisms nb-tensor-prisms
            #self.getnumprisms(),0,0,
            0,0,0,
    #nb-scalar-pyramids nb-vector-pyramids nb-tensor-pyramids
            0,0,0,
    #nb-text2d nb-text2d-chars nb-text3d nb-text3d-chars
            0,0,0,0
            ))
        for timestep in range(len(scaltime)):
            #time-step-values
            f.write("%f "%(timestep*dt))
        f.write("\n")
        #< scalar-triangles-value > ...
        for el in self.elements:
            if el[1] in [pmdtriangle,pmdtrianglerot]: 
                n=(self.getxyz(el[2]),
                self.getxyz(el[3]),
                self.getxyz(el[4]))
                T=[]
                for timestep in range(len(scaltime)):
                    self.scalars=scaltime[timestep]
                    if len(self.nodes) != len(self.scalars):
                        error("Different number of nodes and scalars!")
                    T.extend((self.getscalar(el[2]),
                        self.getscalar(el[3]),
                        self.getscalar(el[4])))
                f.write(formatpos(n,T))
        f.write("$EndView\n")
        f.close()
    def writescalars(self,filename,scalars,C=0.0):
        f=file(filename,"w")
        for n,s in zip(self.nodes,scalars):
            f.write("%4d  %f\n"%(n[0],s+C))
    def writestresspos(self,filename,infotext="PMD_stress"):
        """Writes scalars together with nodes and elements to filename.

        Optional parameter infotext specifies the name of the view 
        in the pos file.

        1) Associates a scalar with every node, so gmsh
        shows points.
        2) Associates a scalar with every node of all elements, so gmsh
        fills the whole element (triangle, quadrangle etc.) with an
        extrapolated color.

        You can set visibility in gmsh (only points, only triangles...).
        """
        if len(self.elements) != len(self.elementdata):
            error("Different number of elements and data!")
        f=file(filename,"w")
        f.write("$PostFormat\n")
        #1.3 file-type data-size
        f.write("%g %d %d\n"%(1.2,0,8))
        f.write("$EndPostFormat\n")
        f.write("$View\n")
        f.write("%s %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d\n"
            "%d %d %d %d\n"%(
    #view-name nb-time-steps
            infotext, 1, 
    #nb-scalar-points nb-vector-points nb-tensor-points
            0,0,0, 
    #nb-scalar-lines nb-vector-lines nb-tensor-lines
            0,0,0,
    #nb-scalar-triangles nb-vector-triangles nb-tensor-triangles
            self.getnumtriangles(),0,0,
    #nb-scalar-quadrangles nb-vector-quadrangles nb-tensor-quadrangles
            self.getnumquadrangles(),0,0,
    #nb-scalar-tetrahedra nb-vector-tetrahedra nb-tensor-tetrahedra
            self.getnumtetrahedra(),0,0,
    #nb-scalar-hexahedra nb-vector-hexahedra nb-tensor-hexahedra
            0,0,0,
    #nb-scalar-prisms nb-vector-prisms nb-tensor-prisms
            0,0,0,
    #nb-scalar-pyramids nb-vector-pyramids nb-tensor-pyramids
            0,0,0,
    #nb-text2d nb-text2d-chars nb-text3d nb-text3d-chars
            0,0,0,0
            ))
        #time-step-values
        f.write("%d\n"%(0))
        #< scalar-triangles-value > ...
        for el,data in zip(self.elements,self.elementdata):
            if el[1] == pmdtriangle: 
                n=(self.getxyz(el[2]),
                self.getxyz(el[3]),
                self.getxyz(el[4]))
                T=data[:3]
                f.write(formatpos(n,T))
        #< scalar-quadrangles-value > ...
        for el,data in zip(self.elements,self.elementdata):
            if el[1] == pmdquadrangle: 
                n=(self.getxyz(el[2]),
                self.getxyz(el[3]),
                self.getxyz(el[4]),
                self.getxyz(el[5]))
                T=data[:4]
                f.write(formatpos(n,T))
        #< scalar-tetrahedra-value > ...
        for el,data in zip(self.elements,self.elementdata):
            if el[1] == pmdtetrahedron: 
                n=(self.getxyz(el[2]),
                self.getxyz(el[3]),
                self.getxyz(el[4]),
                self.getxyz(el[5]))
                T=data[:4]
                f.write(formatpos(n,T))
        f.write("$EndView\n")
        f.close()
    def average(self,list):
        a=0
        for i in list:
            a+=i
        a/=len(list)
        return a
    def average_vectors(self,list):
        a=[0,0,0]
        for i in list:
            a[0]+=i[0]
            a[1]+=i[1]
            a[1]+=i[1]
        a[0]/=len(list)
        a[1]/=len(list)
        a[2]/=len(list)
        return a
    def convert_el_to_nodes(self,els):
        tmp=[]
        for i in self.nodes:
            tmp.append([])
        for el,data in zip(self.elements,els):
            nodes=el[2:]
            if self.is2d:
                assert(data[len(nodes)*2]==0)
            for node in nodes:
                tmp[node-1].append(data)  
        self.scalars=[]
        for s in tmp:
            assert s!=[]
            self.scalars.append(self.average(s))
    def convert_stress_to_nodes(self):
        tmp=[]
        for i in self.nodes:
            tmp.append([])
        for el,data in zip(self.elements,self.elementdata):
            nodes=el[2:]
            if self.is2d:
                assert(data[len(nodes)*2]==0)
            for node,scalar in zip(nodes,data):
                tmp[node-1].append(scalar)  
        self.scalars=[]
        for s in tmp:
            self.scalars.append(self.average(s))
    def getxyz(self,n):
        "Returns a tuple (x,y,z) of a node whose number is n."
        n-=1
        if n<0 or n>=len(self.nodes):
            error("node not found")
        node=self.nodes[n]
        if node[0] != n+1: error("node numbers mischmatch",3)
        return [node[1],node[2],node[3]]
    def getscalar(self,n):
        "Returns a scalar of a node whose number is n."
        n-=1
        if n<0 or n>=len(self.scalars):
            error("node not found")
        return self.scalars[n]
    def getvector(self,n):
        "Returns a vector associated to a node whose number is n."
        n-=1
        if n<0 or n>=len(self.vectors):
            error("node not found")
        return self.vectors[n]
    def getnumtriangles(self):
        n=0
        for x in self.elements: 
            if x[1] in [pmdtriangle,pmdtrianglerot]: n+=1
        return n
    def getnumquadrangles(self):
        n=0
        for x in self.elements: 
            if x[1] == pmdquadrangle: n+=1
        return n
    def getnumtetrahedra(self):
        n=0
        for x in self.elements: 
            if x[1] == pmdtetrahedron: n+=1
        return n
    def getnumhexahedra(self):
        n=0
        for x in self.elements: 
            if x[1] == pmdhexahedron: n+=1
        return n
    def getnumprisms(self):
        n=0
        for x in self.elements: 
            if x[1] == pmdprism: n+=1
        return n
    def removecentralnodes(self,nods):
        if nods==[]:
            return
        nods.sort()
        n=nods[0]
        for nod in nods:
            if n!=nod:
                #n=nod
                error("faces aren't consecutive",2)
            n+=1
        if nods[-1] != self.nodes[-1][0]:
            error("faces aren't at the end of self.nodes",2)
        self.nodes=self.nodes[:len(self.nodes)-len(nods)]
    def getscalar2(self,n,scalars):
        "Returns a scalar of a node whose number is n."
        n-=1
        if n<0 or n>=len(scalars):
            error("node not found")
        return scalars[n]
    def scalars_elements2nodes(self,scalarsel):
        assert(len(self.elements)==len(scalarsel))
        tmp=[]
        for i in self.nodes:
            tmp.append([])
        for el,data in zip(self.elements,scalarsel):
            nodes=el[2:]
            #if self.is2d:
            #   assert(data[len(nodes)*2]==0)
            for node,scalar in zip(nodes,data):
                tmp[node-1].append(scalar)  
        scalars=[]
        for s in tmp:
            scalars.append(self.average(s))
        return scalars
    def vectors_elements2nodes(self,scalarsel):
        assert(len(self.elements)==len(scalarsel))
        tmp=[]
        for i in self.nodes:
            tmp.append([])
        for el,data in zip(self.elements,scalarsel):
            nodes=el[2:]
            #if self.is2d:
            #   assert(data[len(nodes)*2]==0)
            for node,scalar in zip(nodes,data):
                tmp[node-1].append(scalar)  
        scalars=[]
        for s in tmp:
            if s==[]: error("There are some extra nodes!")
            scalars.append(self.average_vectors(s))
        return scalars
    def dist(self,a,b):
        p=self.getxyz(a)
        q=self.getxyz(b)
        return math.sqrt((p[0]-q[0])**2+(p[1]-q[1])**2+(p[2]-q[2])**2)
    def dist2(self,p,q):
        return math.sqrt((p[0]-q[0])**2+(p[1]-q[1])**2+(p[2]-q[2])**2)
    def det(self,x,y):
        return x[0]*y[1]+x[1]*y[2]+x[2]*y[0]-(x[0]*y[2]+x[1]*y[0]+x[2]*y[1])

    def computegrad(self,scalars):
        """Returns the gradient of ``scalars`` (both are given in nodes)."""
        grad=[]
        for e in self.elements:
            if e[1]==pmdtriangle or e[1]==pmdtrianglerot:
                #nodes e[2],e[3],e[4]
                a1=self.getxyz(e[2])
                a1[2]=self.getscalar2(e[2],scalars)
                a2=self.getxyz(e[3])
                a2[2]=self.getscalar2(e[3],scalars)
                a3=self.getxyz(e[4])
                a3[2]=self.getscalar2(e[4],scalars)

                a=self.det((a1[1],a2[1],a3[1]),(a1[2],a2[2],a3[2]))
                b=-self.det((a1[0],a2[0],a3[0]),(a1[2],a2[2],a3[2]))
                c=self.det((a1[0],a2[0],a3[0]),(a1[1],a2[1],a3[1]))
                v=(-a/c,-b/c,0)
                #print v
                grad.append((v,v,v))
            elif e[1]==pmdtetrahedron:
                #quick hack...
                #nodes e[2],e[3],e[4],e[5]
                a1=self.getxyz(e[2])
                a2=self.getxyz(e[3])
                d=(self.getscalar2(e[2],scalars)-self.getscalar2(e[3],scalars))\
                    /self.dist2(a1,a2)
                v=((a1[0]-a2[0])*d,
                    (a1[1]-a2[1])*d,
                    (a1[2]-a2[2])*d)
                #if e[0]==2879: print v
                grad.append((v,v,v,v))
            else:
                error("Element not implemented yet.")
        return self.vectors_elements2nodes(grad)

    def computenorm(self,vectors):
        """Returns sqrt(x^2+y^2+z^2) for all vectors (x,y,z) in ``vectors``."""
        scalars=[]
        for v in vectors:
            scalars.append(self.dist2(v,(0,0,0)))
        return scalars
            
    def printinfo(self):
        print "nodes:    %d"%(len(self.nodes))
        print "elements: %d"%(len(self.elements))

    def readGMV(self,filename,what=2):
        """Reads GMV file.

        what ... 0 read only mesh
            ... 1 read only data
            ... 2 read both
        """
        if what in [0,2]:
            self.clean()
        f=file(filename)
        l=f.readline(); check(l,"gmvinput ascii\n")
        l=f.readline(); check(l,"\n")
        l=f.readline(); 
        x=l.split()
        check(x[0],"nodes")
        nnod=int(x[1])
        l=f.readline(); 
        if what in [0,2]: 
            xs=[float(x) for x in l.split()]
            check(len(xs),nnod)
        l=f.readline(); 
        if what in [0,2]: 
            ys=[float(x) for x in l.split()]
            check(len(ys),nnod)
        l=f.readline(); 
        if what in [0,2]: 
            zs=[float(x) for x in l.split()]
            check(len(zs),nnod)
            nodes=zip(range(1,nnod+1),xs,ys,zs)
        l=f.readline(); 
        check(l,"\n")
        l=f.readline(); 
        x=l.split()
        check(x[0],"cells")
        if what in [0,2]:
            nelm=int(x[1])
            elements=[]
            for i in range(1,nelm+1):
                l=f.readline(); 
                x=l.split()
                if x[0]=="quad":
                    check(x[1],"4")
                    l=f.readline(); 
                    p=[int(x) for x in l.split()]
                    elements.append((i,pmdquadrangle)+tuple(p))
                elif x[0]=="phex8":
                    check(x[1],"8")
                    l=f.readline(); 
                    p=[int(x) for x in l.split()]
                    elements.append((i,pmdhexahedron)+tuple(p))
                elif x[0]=="tri":
                    check(x[1],"3")
                    l=f.readline(); 
                    p=[int(x) for x in l.split()]
                    elements.append((i,pmdtriangle)+tuple(p))
                elif x[0]=="tet":
                    check(x[1],"4")
                    l=f.readline(); 
                    p=[int(x) for x in l.split()]
                    elements.append((i,pmdtetrahedron)+tuple(p))
                else:
                    error("unimplemented yet.",2)
            self.nodes=nodes
            self.elements=elements
        if what in [1,2]:
            while l!="variable\n":
                l=f.readline(); 
            l=f.readline(); 
            l=f.readline(); 
            scalars=[float(x) for x in l.split()]
            check(len(scalars),nnod)
            return scalars
    def writeregions(self,filename):
        f=open(filename,"w")
        f.write(str(self.regions))
    def writeBC(self,filename,verbose=True):
        """self.faces contain triplets (p1,p2,p3) which are triangles of
        tetrahedrons on the boundary. We need to find the number of each
        corresponding tetrahedron and it's side."""
        def findelements(faces,elements):
            """Returns element numbers (in a list), which lies
            at the boundary determined by faces.
            
            Also returns the side of element. Currently only support
            triangles and quadrangles.
            """
            t=len(elements)
            c=0
            el=[]
            for p in elements:
                c+=1
                if c%500==0: print 100.0*c/t
                for ii,f in enumerate(faces):
                    nods=[]
                    for i,n in enumerate(p[2:]): 
                        if n in f: nods.append(i+1)
                    if len(nods)!=len(f): continue
                    if len(f)==3:
                        if p[1] == pmdtetrahedron:
                            if nods==[1,2,3]:
                                side=1
                            elif nods==[1,2,4]:
                                side=2
                            elif nods==[2,3,4]:
                                side=3
                            elif nods==[1,3,4]:
                                side=4
                            else:
                                raise"findelements: tetrahedron face mischmatch"
                            el.append((p[0],side))
                            del faces[ii]
                            break;
                        else:
                            raise "findelements: unsupported element in mesh"
                    else:
                        raise "findelements: unsupported face %s"%(repr(f))
            return el
        def buildmapping(elements):
            m={}
            for e in elements:
                for n in e[2:]:
                    if not m.has_key(n): m[n]=[]
                    m[n].append(e[0])
            return m
        def findelements2(faces,elements,nemap):
            bc=[]
            for f in faces:
                assert len(f)==3
                candidates=set(nemap[f[0]])
                candidates.intersection_update(nemap[f[1]])
                candidates.intersection_update(nemap[f[2]])
                assert len(candidates)==1  #the face "f" belongs to just 1 el.
                elnum=candidates.pop()
                el=elements[elnum-1]
                assert el[0]==elnum  #the mapping "nemap" is correct
                nods=[]
                for i,n in enumerate(el[2:]): 
                    if n in f: nods.append(i+1)
                assert len(nods)==len(f) #the mapping "nemap" is correct
                if el[1] == pmdtetrahedron:
                    if nods==[1,2,3]:
                        side=1
                    elif nods==[1,2,4]:
                        side=2
                    elif nods==[2,3,4]:
                        side=3
                    elif nods==[1,3,4]:
                        side=4
                    else:
                        raise"findelements: tetrahedron face mischmatch"
                    bc.append((elnum,side))
                else:
                    raise "findelements: unsupported element in mesh"
            return bc

        up=progressbar.MyBar("Writing BC to %s:"%filename)
        if verbose: up.init(2*len(self.faces))
        nemap=buildmapping(self.elements)
        bc={}
        c=0
        for key in self.faces:
            #bc[key]=findelements(self.faces[key],self.elements)
            bc[key]=findelements2(self.faces[key],self.elements,nemap)
            c+=1
            if verbose: up.update(c)
        f=open(filename,"w")
        #f.write(repr(bc))
        f.write("%d\n"%len(bc))
        for k in bc:
            f.write("%d %d %s\n"%(k,len(bc[k]),numlist2str(flat(bc[k]))))
            c+=1
            if verbose: up.update(c)
#        print bc[2]
#        for i in range(len(bc[2])):
#            print bc[2][i], self.elements[bc[2][i][0]-1] 

def flat(a):
    r=[]
    for x in a:
        if isinstance(x,list) or isinstance(x,tuple):
            r.extend(flat(x))
        else:
            r.append(x)
    return r

def numlist2str(x):
    s=""
    for i in x:
        s+="%d "%i
    return s[:-1]

def formatpos(n,T):
    #n=(n1,n2,n3...)
    #T=(T1,T2,T3...)
    #nodes of the triangle (for example):
    #n1=(x1,y1,z1), n2=(x2,y2,z2) and n3=(x3,y3,z3)
    #format of *.pos:
    #x1,x2,x3,y1,y2,y3,z1,z2,z3,T1,T2,T3 
    y=[]
    for i in range(len(n[0])):
        for j in range(len(n)):
            y.append(n[j][i])
    y.extend(T)
    str="%.18f "*(len(y))%tuple(y)
    return "%s\n"%(str)
def formatpos2(n,T):
    #n=(n1,n2,n3...)
    #T=(T1,T2,T3...)
    #nodes of the triangle (for example):
    #n1=(x1,y1,z1), n2=(x2,y2,z2) and n3=(x3,y3,z3)
    #format of *.pos:
    #x1,x2,x3,y1,y2,y3,z1,z2,z3,T1,T2,T3 
    y=[]
    for i in range(len(n[0])):
        for j in range(len(n)):
            y.append(n[j][i])
    for v in T:
        y.extend(v)
    str="%.18f "*(len(y))%tuple(y)
    return "%s\n"%(str)

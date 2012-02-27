from sfepy.fem import Mesh
import numpy as nm

def meshgen_from_goem(geo, a=None, quadratic=False, verbose=True,
                      refine=False, polyfilename='./meshgen.poly',
                      out='mesh', **kwargs):
    """
    Runs mesh generator - tetgen for 3D or triangle for 2D meshes.

    Parameters
    ----------
    geo : geometry
        geometry description
    a : int, optional
        a maximum area/volume constraint
    quadratic : bool, optional
        set True for quadratic elements
    verbose : bool, optional
        detailed information
    refine : bool, optional
        refines mesh

    Returns
    -------
    mesh : mesh
        triangular or tetrahedral mesh
    """

    import os.path as op
    import pexpect

    # write geometry to poly file
    geo.to_poly_file(polyfilename)

    if not refine:
        params = "-Apq"
    else:
        params = "-Arq"
    if verbose:
        params = params + " -Q"
    if a != None and not refine:
        params = params + " -a%f" % (a)
    if refine:
        params = params + " -a"
    if quadratic:
        params = params + " -o2"
    params = params + " %s" % (polyfilename)

    meshgen_call = {2: 'triangle', 3: 'tetgen'}
    cmd = "%s %s" % (meshgen_call[geo.dim], params)
    if verbose: print "Generating mesh using", cmd
    if geo.dim == 2:
        p=pexpect.run(cmd, timeout=None)
        bname, ext = op.splitext(polyfilename)
        mesh = Mesh.from_file(bname + '.1.node')
        mesh.write(bname + '.' + out)
    if geo.dim == 3:
        p=pexpect.spawn(cmd, timeout=None)
        if not refine:
            p.expect("Opening %s." % (polyfilename))
        else:
            p.expect("Opening %s.node.\r\n" % (polyfilename))
            p.expect("Opening %s.ele.\r\n" % (polyfilename))
            p.expect("Opening %s.face.\r\n" % (polyfilename))
            p.expect("Opening %s.vol." % (polyfilename))
        assert p.before == ""
        p.expect(pexpect.EOF)
        if p.before != "\r\n":
            print p.before
            raise "Error when running mesh generator (see above for output): %s" % cmd

# http://www.cs.cmu.edu/~quake/triangle.html
#
# triangle [-prq__a__uAcDjevngBPNEIOXzo_YS__iFlsCQVh] input_file
#     -p  Triangulates a Planar Straight Line Graph (.poly file).
#     -r  Refines a previously generated mesh.
#     -q  Quality mesh generation.  A minimum angle may be specified.
#     -a  Applies a maximum triangle area constraint.
#     -u  Applies a user-defined triangle constraint.
#     -A  Applies attributes to identify triangles in certain regions.
#     -c  Encloses the convex hull with segments.
#     -D  Conforming Delaunay:  all triangles are truly Delaunay.
#     -j  Jettison unused vertices from output .node file.
#     -e  Generates an edge list.
#     -v  Generates a Voronoi diagram.
#     -n  Generates a list of triangle neighbors.
#     -g  Generates an .off file for Geomview.
#     -B  Suppresses output of boundary information.
#     -P  Suppresses output of .poly file.
#     -N  Suppresses output of .node file.
#     -E  Suppresses output of .ele file.
#     -I  Suppresses mesh iteration numbers.
#     -O  Ignores holes in .poly file.
#     -X  Suppresses use of exact arithmetic.
#     -z  Numbers all items starting from zero (rather than one).
#     -o2 Generates second-order subparametric elements.
#     -Y  Suppresses boundary segment splitting.
#     -S  Specifies maximum number of added Steiner points.
#     -i  Uses incremental method, rather than divide-and-conquer.
#     -F  Uses Fortune's sweepline algorithm, rather than d-and-c.
#     -l  Uses vertical cuts only, rather than alternating cuts.
#     -s  Force segments into mesh by splitting (instead of using CDT).
#     -C  Check consistency of final mesh.
#     -Q  Quiet:  No terminal output except errors.
#     -V  Verbose:  Detailed information on what I'm doing.
#     -h  Help:  Detailed instructions for Triangle.

# http://tetgen.berlios.de/
#
# tetgen [-prq_a_AiMYS_T_dzo_fenvgGOJBNEFICQVh] input_file
#     -p  Tetrahedralizes a piecewise linear complex (PLC).
#     -r  Reconstructs a previously generated mesh.
#     -q  Refines mesh (to improve mesh quality).
#     -a  Applies a maximum tetrahedron volume constraint.
#     -A  Assigns attributes to tetrahedra in different regions.
#     -i  Inserts a list of additional points into mesh.
#     -M  No merge of coplanar facets.
#     -Y  No splitting of input boundaries (facets and segments).
#     -S  Specifies maximum number of added points.
#     -T  Sets a tolerance for coplanar test (default 1e-8).
#     -d  Detects self-intersections of facets of the PLC.
#     -z  Numbers all output items starting from zero.
#     -o2 Generates second-order subparametric elements.
#     -f  Outputs all faces to .face file.
#     -e  Outputs all edges to .edge file.
#     -n  Outputs tetrahedra neighbors to .neigh file.
#     -v  Outputs Voronoi diagram to files.
#     -g  Outputs mesh to .mesh file for viewing by Medit.
#     -G  Outputs mesh to .msh file for viewing by Gid.
#     -O  Outputs mesh to .off file for viewing by Geomview.
#     -K  Outputs mesh to .vtk file for viewing by Paraview.
#     -J  No jettison of unused vertices from output .node file.
#     -B  Suppresses output of boundary information.
#     -N  Suppresses output of .node file.
#     -E  Suppresses output of .ele file.
#     -F  Suppresses output of .face file.
#     -I  Suppresses mesh iteration numbers.
#     -C  Checks the consistency of the final mesh.
#     -Q  Quiet:  No terminal output except errors.
#     -V  Verbose:  Detailed information, more terminal output.
#     -h  Help:  A brief instruction for using TetGen.

def elems_q2t(el):

    nel, nnd = el.shape
    if nnd > 4:
        q2t = nm.array([[0, 2, 3, 6],
                        [0, 3, 7, 6],
                        [0, 7, 4, 6],
                        [0, 5, 6, 4],
                        [1, 5, 6, 0],
                        [1, 6, 2, 0]])

    else:
        q2t = nm.array([[0, 1, 2],
                        [0, 2, 3]])

    ns, nn = q2t.shape
    nel *= ns

    out = nm.zeros((nel, nn), dtype=nm.int32);

    for ii in range(ns):
        idxs = nm.arange(ii, nel, ns)

        out[idxs,:] = el[:, q2t[ii,:]]

    return nm.ascontiguousarray(out)

def meshgen_from_voxels(voxels, dims, etype='q'):
    """
    Generate FE mesh from voxels (volumetric data).

    Parameters
    ----------
    voxels : array
        Voxel matrix, 1=material.
    dims : array
        Size of one voxel.
    etype : integer, optional
        'q' - quadrilateral or hexahedral elements
        't' - triangular or tetrahedral elements
    Returns
    -------
    mesh : mesh
        Finite element mesh.
    """

    dims = dims.squeeze()
    dim = len(dims)
    nddims = nm.array(voxels.shape) + 1

    nodemtx = nm.zeros(nddims, dtype=nm.int32)

    if dim == 2:
        #iy, ix = nm.where(voxels.transpose())
        iy, ix = nm.where(voxels)
        nel = ix.shape[0]

        if etype == 'q':
            nodemtx[ix,iy] += 1
            nodemtx[ix + 1,iy] += 1
            nodemtx[ix + 1,iy + 1] += 1
            nodemtx[ix,iy + 1] += 1

        elif etype == 't':
            nodemtx[ix,iy] += 2
            nodemtx[ix + 1,iy] += 1
            nodemtx[ix + 1,iy + 1] += 2
            nodemtx[ix,iy + 1] += 1
            nel *= 2

    elif dim == 3:
        #iy, ix, iz = nm.where(voxels.transpose(1, 0, 2))
        iy, ix, iz = nm.where(voxels)
        nel = ix.shape[0]

        if etype == 'q':
            nodemtx[ix,iy,iz] += 1
            nodemtx[ix + 1,iy,iz] += 1
            nodemtx[ix + 1,iy + 1,iz] += 1
            nodemtx[ix,iy + 1,iz] += 1
            nodemtx[ix,iy,iz + 1] += 1
            nodemtx[ix + 1,iy,iz + 1] += 1
            nodemtx[ix + 1,iy + 1,iz + 1] += 1
            nodemtx[ix,iy + 1,iz + 1] += 1

        elif etype == 't':
            nodemtx[ix,iy,iz] += 6
            nodemtx[ix + 1,iy,iz] += 2
            nodemtx[ix + 1,iy + 1,iz] += 2
            nodemtx[ix,iy + 1,iz] += 2
            nodemtx[ix,iy,iz + 1] += 2
            nodemtx[ix + 1,iy,iz + 1] += 2
            nodemtx[ix + 1,iy + 1,iz + 1] += 6
            nodemtx[ix,iy + 1,iz + 1] += 2
            nel *= 6

    else:
        msg = 'incorrect voxel dimension! (%d)' % dim
        raise ValueError(msg)

    ndidx = nm.where(nodemtx)
    coors = nm.array(ndidx).transpose() * dims
    nnod = coors.shape[0]

    nodeid = -nm.ones(nddims, dtype=nm.int32)
    nodeid[ndidx] = nm.arange(nnod)

    # generate elements
    if dim == 2:
        elems = nm.array([nodeid[ix,iy],
                          nodeid[ix + 1,iy],
                          nodeid[ix + 1,iy + 1],
                          nodeid[ix,iy + 1]]).transpose()

    elif dim == 3:
        elems = nm.array([nodeid[ix,iy,iz],
                          nodeid[ix + 1,iy,iz],
                          nodeid[ix + 1,iy + 1,iz],
                          nodeid[ix,iy + 1,iz],
                          nodeid[ix,iy,iz + 1],
                          nodeid[ix + 1,iy,iz + 1],
                          nodeid[ix + 1,iy + 1,iz + 1],
                          nodeid[ix,iy + 1,iz + 1]]).transpose()

    if etype == 't':
        elems = elems_q2t(elems)

    eid = etype + str(dim)
    eltab = {'q2': 4, 'q3': 8, 't2': 3, 't3': 4}

    mesh = Mesh.from_data('voxel_data',
                          coors, nm.ones((nnod, ), dtype=nm.int32),
                          {0: nm.ascontiguousarray(elems)},
                          {0: nm.ones((nel, ), dtype=nm.int32)},
                          {0: '%d_%d' % (dim, eltab[eid])})

    return mesh

def meshgen_from_poly(filename, verbose=True):
    """
    Import mesh generated by tetgen or triangle.

    Parameters
    ----------
    filename : string
        file name

    Returns
    -------
    mesh : mesh
        triangular or tetrahedral mesh
    """

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
    m=Mesh()
    m.nodes=getnodes(filename+".node")
    m.elements,m.regions, lin=getele(filename+".ele")
    if not lin:
        #tetgen doesn't compute xyz coordinates of the aditional 6 nodes
        #(only of the 4 corner nodes) in tetrahedra.
        calculatexyz(m.nodes,m.elements)
    m.faces=getBCfaces(filename+".face")
    return m


def smooth_mesh(mesh, n_iter=4, lam=0.6307, mu=-0.6347,
                weights=None, bconstr=True,
                volume_corr=False):
    """
    FE mesh smoothing.

    Based on:

    [1] Steven K. Boyd, Ralph Muller, Smooth surface meshing for automated
    finite element model generation from 3D image data, Journal of
    Biomechanics, Volume 39, Issue 7, 2006, Pages 1287-1295,
    ISSN 0021-9290, 10.1016/j.jbiomech.2005.03.006.
    (http://www.sciencedirect.com/science/article/pii/S0021929005001442)

    Parameters
    ----------
    mesh : mesh
        FE mesh.
    n_iter : integer, optional
        Number of iteration steps.
    lam : float, optional
        Smoothing factor, see [1].
    mu : float, optional
        Unshrinking factor, see [1].
    weights : array, optional
        Edge weights, see [1].
    bconstr: logical, optional
        Boundary constraints, if True only surface smoothing performed.
    volume_corr: logical, optional
        Correct volume after smoothing process.

    Returns
    -------
    coors : array
        Coordinates of mesh nodes.
    """

    from sfepy.fem import Domain
    import scipy.sparse as sps
    import scipy as sc

    def laplacian(coors, weights):

        n_nod = coors.shape[0]
        displ = (weights - sps.identity(n_nod)) * coors

        return displ

    def taubin(coors0, weights, lam, mu, n_iter):

        coors = coors0.copy()

        for ii in range(n_iter):
            displ = laplacian(coors, weights)
            if nm.mod(ii, 2) == 0:
                coors += lam * displ
            else:
                coors += mu * displ

        return coors

    def get_volume(el, nd):

        dim = nd.shape[1]
        nnd = el.shape[1]

        etype = '%d_%d' % (dim, nnd)
        if etype == '2_4' or etype == '3_8':
            el = elems_q2t(el)

        vol = 0.0
        bc = nm.zeros((dim, ), dtype=nm.double)
        mtx = nm.ones((dim + 1, dim + 1), dtype=nm.double)
        mul = 1.0 / sc.factorial(dim)
        if dim == 3:
            mul *= -1.0

        for iel in el:
            mtx[:,:-1] = nd[iel,:]
            ve = mul * nm.linalg.det(mtx)
            vol += ve
            bc += ve * mtx.sum(0)[:-1] / nnd

        bc /= vol

        return vol, bc


    domain = Domain('mesh', mesh)

    n_nod = mesh.n_nod
    edges = domain.ed

    if weights is None:
        # initiate all vertices as inner - hierarchy = 2
        node_group = nm.ones((n_nod,), dtype=nm.int16) * 2
        # boundary vertices - set hierarchy = 4
        if bconstr:
            # get "nodes of surface"
            if domain.fa: # 3D.
                fa = domain.fa
            else:
                fa = domain.ed

            flag = fa.mark_surface_facets()
            ii = nm.where( flag > 0 )[0]
            aux = nm.unique(fa.facets[ii])
            if aux[0] == -1: # Triangular faces have -1 as 4. point.
                aux = aux[1:]

            node_group[aux] = 4

        # generate costs matrix
        costs = sps.lil_matrix((n_nod, n_nod), dtype=nm.double)
        for ied in range(edges.mtx.shape[0]):
            cc = edges.mtx.getrow(ied).indices[0]
            n1, n2 = edges.facets[cc,:]
            if node_group[n2] >= node_group[n1]:
                costs[n1, n2] = 1.0

            if node_group[n1] >= node_group[n2]:
                costs[n2, n1] = 1.0

        # generate weights matrix
        aux = sps.lil_matrix((n_nod, n_nod), dtype=nm.double)
        aux.setdiag(1.0 / costs.sum(1))
        weights = (aux * costs).tocsr()

    coors = taubin(mesh.coors, weights, lam, mu, n_iter)

    if volume_corr:

        volume0, bc = get_volume(mesh.conns[0], mesh.coors)
        volume, _ = get_volume(mesh.conns[0], coors)

        scale = volume0 / volume

        coors = (coors - bc) * scale + bc

    return coors

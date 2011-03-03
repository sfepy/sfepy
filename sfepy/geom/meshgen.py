import pexpect
from meshgenutils import write_poly
from sfepy.fem import Mesh
import os.path as op

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

    # write geometry to poly file
    write_poly(geo, polyfilename)

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

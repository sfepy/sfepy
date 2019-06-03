import numpy as np
import sys
import gmsh


def triangle_max_edge(x):
    a = np.sum((x[:,0,:]-x[:,1,:])**2,1)**0.5
    b = np.sum((x[:,0,:]-x[:,2,:])**2,1)**0.5
    c = np.sum((x[:,1,:]-x[:,2,:])**2,1)**0.5
    return np.maximum(a,np.maximum(b,c))


class Mesh:

    def __init__(self):
        self.vtags, vxyz, _ = gmsh.model.mesh.getNodes()
        self.vxyz = vxyz.reshape((-1,3))
        vmap = dict({j:i for i,j in enumerate(self.vtags)})
        self.triangles_tags, evtags = gmsh.model.mesh.getElementsByType(2)
        evid = np.array([vmap[j] for j in evtags])
        self.triangles = evid.reshape((self.triangles_tags.shape[-1],-1))


def my_function(xyz):
    a = 6*(np.hypot(xyz[...,0]-.5,xyz[...,1]-.5)-.2)
    f = np.real(np.arctanh(a+0j))
    return f


def compute_interpolation_error(nodes, triangles, f):
    jac,det,pt = gmsh.model.mesh.getJacobians(2,"Gauss2")
    uvwo,numcomp,sf = gmsh.model.mesh.getBasisFunctions(2,"Gauss2","Lagrange")
    weights = uvwo.reshape([-1,4])[:,3]
    sf = sf.reshape((weights.shape[0],-1))
    qx = pt.reshape((triangles.shape[0],-1,3))
    det = np.abs(det.reshape((triangles.shape[0],-1)))
    f_vert = f(nodes)
    f_fem = np.dot(f_vert[triangles],sf)
    err_tri = np.sum((f_fem-f(qx))**2*det*weights,1)
    return f_vert, np.sqrt(err_tri)


def compute_size_field(nodes, triangles, err, N):
    x = nodes[triangles]
    a = 2.
    d = 2.
    fact = (a**((2.+a)/(1.+a)) + a**(1./(1.+a))) * np.sum(err**(2./(1.+a)))
    ri = err**(2./(2.*(1+a))) * a**(1./(d*(1.+a))) * ((1.+a)*N/fact)**(1./d)
    return triangle_max_edge(x)/ri


print("Usage: adapt_mesh [intial lc] [target #elements] [dump files]")

lc = 0.02;
N = 10000;
dumpfiles = False
if len(sys.argv) > 1: lc = float(sys.argv[1])
if len(sys.argv) > 2: N = int(sys.argv[2])
if len(sys.argv) > 3: dumpfiles = int(sys.argv[3])

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

# create a geometrical gmsh.model
gmsh.model.add("square")
square = gmsh.model.occ.addRectangle(0, 0, 0, 1, 1)
gmsh.model.occ.synchronize()

# create intial uniform mesh
pnts = gmsh.model.getBoundary([(2,square)], True, True, True);
gmsh.model.mesh.setSize(pnts, lc)
gmsh.model.mesh.generate(2)
if dumpfiles: gmsh.write("mesh.msh")
mesh = Mesh()

# compute and visualize the interpolation error
f_nod, err_ele = compute_interpolation_error(mesh.vxyz, mesh.triangles, my_function)
f_view = gmsh.view.add("nodal function")
gmsh.view.addModelData(f_view, 0, "square", "NodeData",
                       mesh.vtags, f_nod[:,None])
if dumpfiles: gmsh.view.write(f_view, "f.pos")
err_view = gmsh.view.add("element-wise error")
gmsh.view.addModelData(err_view, 0, "square", "ElementData",
                       mesh.triangles_tags, err_ele[:,None])
if dumpfiles: gmsh.view.write(err_view, "err.pos")

# compute and visualize the remeshing size field
sf_ele = compute_size_field(mesh.vxyz,mesh.triangles, err_ele, N)
sf_view = gmsh.view.add("mesh size field")
gmsh.view.addModelData(sf_view, 0, "square", "ElementData",
                       mesh.triangles_tags, sf_ele[:,None])
if dumpfiles: gmsh.view.write(sf_view, "sf.pos")

# create a new gmsh.model (to remesh the original gmsh.model in-place, the size field
# should be created as a list-based view)
gmsh.model.add("square2")
gmsh.model.occ.addRectangle(0, 0, 0, 1, 1)
gmsh.model.occ.synchronize()

# mesh the new gmsh.model using the size field
bg_field = gmsh.model.mesh.field.add("PostView");
gmsh.model.mesh.field.setNumber(bg_field, "ViewTag", sf_view);
gmsh.model.mesh.field.setAsBackgroundMesh(bg_field);
gmsh.model.mesh.generate(2)
if dumpfiles: gmsh.write("mesh2.msh")
mesh2 = Mesh()

# compute and visualize the interpolation error on the adapted mesh
f2_nod, err2_ele = compute_interpolation_error(mesh2.vxyz,mesh2.triangles, my_function)
f2_view = gmsh.view.add("nodal function on adapted mesh")
gmsh.view.addModelData(f2_view, 0, "square2", "NodeData",
                       mesh2.vtags, f2_nod[:,None])
if dumpfiles: gmsh.view.write(f2_view, "f2.pos")
err2_view = gmsh.view.add("element-wise error on adapated mesh")
gmsh.view.addModelData(err2_view, 0, "square2", "ElementData",
                       mesh2.triangles_tags, err2_ele[:,None])
if dumpfiles: gmsh.view.write(err2_view, "err2.pos")

# show everything in the gui
gmsh.fltk.run()

gmsh.finalize()

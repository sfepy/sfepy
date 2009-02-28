from sfepy.base.base import *
from sfepy.fem.mesh import make_inverse_connectivity, find_closest_nodes
from sfepy.base.la import inverse_element_mapping

class Probe(Struct):

    def __call__(self, variable):
        return self.probe(variable)

class LineProbe(Probe):

    def __init__(self, p0, p1, n_point, mesh):
        p0 = nm.array(p0, dtype=nm.float64)
        p1 = nm.array(p1, dtype=nm.float64)
        Probe.__init__(self, p0=p0, p1=p1, n_point=n_point, mesh=mesh)

        iconn = make_inverse_connectivity(mesh.conns, mesh.n_nod,
                                          combine_groups=True)
        self.iconn = iconn

        dirvec = self.p1 - self.p0
        self.length = nm.linalg.norm(dirvec)
        self.dirvec = dirvec / self.length

    def probe(self, variable):
        print variable

        field = variable.field
        vdim = field.dim[0]

        vals = nm.empty((vdim, self.n_point), dtype=variable.dtype)

        coor = self.mesh.nod0[:,:-1]
        conns = self.mesh.conns
        tt = time.clock()
        for ii, step in enumerate(nm.linspace(0, self.length, self.n_point)):
            point = nm.array(self.p0 + self.dirvec * step, ndmin=2)
##             print ii, step, point
            # This could go in Variable.interp_to_point().
            ic = find_closest_nodes(coor, point)
            els = self.iconn[ic]
            for ig, iel in els:
##                 print ig, iel
                nodes = conns[ig][iel]
                ecoor = coor[nodes]
##                 print ecoor

                interp = field.aps[ig].interp
                base_fun = interp.base_funs['v'].fun

                xi = inverse_element_mapping(point, ecoor, base_fun,
                                             suppress_errors=True)
                try:
                    # Verify that we are inside the element.
                    bf = base_fun.value(xi[nm.newaxis,:], base_fun.nodes,
                                        suppress_errors=False)
                except AssertionError:
                    continue
                break

##             print xi, bf

            # For scalar fields only!!!
            vals[:,ii] = nm.dot(bf,variable()[nodes])
        print time.clock() - tt
        print vals
        import pylab
        pylab.plot(vals.squeeze())
        pylab.show()

"""
This example shows how to use the post_process_hook and probe_hook options.

Use it as follows (assumes running from the sfepy directory; on Windows, you
may need to prefix all the commands with "python " and remove "./"):

1. solve the problem:

   ./simple.py examples/linear_elasticity/linear_elastic_probes.py

2. optionally, view the results:

   ./postproc.py cylinder.h5 -b

3. optionally, convert results to VTK, and view again:

   ./extractor.py -d cylinder.h5
   ./postproc.py cylinder.vtk -b

4. probe the data:

   ./probe.py examples/linear_elasticity/linear_elastic_probes.py cylinder.h5

Find :math:`\ul{u}` such that:

.. math::
    \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})
    = 0
    \;, \quad \forall \ul{v} \;,

where

.. math::
    D_{ijkl} = \mu (\delta_{ik} \delta_{jl}+\delta_{il} \delta_{jk}) +
    \lambda \ \delta_{ij} \delta_{kl}
    \;.
"""
# Just grab the problem definition of linear_elastic.py.
from linear_elastic import *

# Define options.
options = {
    'output_dir' : '.',
    'output_format' : 'h5', # VTK reader cannot read cell data yet...

    'post_process_hook' : 'post_process',
    'gen_probes' : 'gen_lines',
    'probe_hook' : 'probe_hook',
}

# Update materials, as ev_cauchy_stress below needs the elastic constants in
# the tensor form.
from sfepy.mechanics.matcoefs import stiffness_from_lame

solid = materials['solid'][0]
lam, mu = solid['lam'], solid['mu']
solid.update({
    'D' : stiffness_from_lame(3, lam=lam, mu=mu),
})

# Update fields and variables to be able to use probes for tensors.
fields.update({
    'sym_tensor': ('real', 6, 'Omega', 0),
})

variables.update({
    's' : ('parameter field', 'sym_tensor', None),
})

# Define the function post_process, that will be called after the problem is
# solved.
def post_process(out, problem, state, extend=False):
    """
    This will be called after the problem is solved.

    Parameters
    ----------
    out : dict
        The output dictionary, where this function will store additional data.
    problem : ProblemDefinition instance
        The current ProblemDefinition instance.
    state : array
        The computed state vector, containing FE coefficients of all the
        unknown variables.
    extend : bool
        The flag indicating whether to extend the output data to the whole
        domain. It can be ignored if the problem is solved on the whole domain
        already.

    Returns
    -------
    out : dict
        The updated output dictionary.
    """
    from sfepy.base.base import Struct

    # Cauchy strain averaged in elements.
    strain = problem.evaluate('ev_cauchy_strain.i1.Omega( u )',
                              mode='el_avg')
    out['cauchy_strain'] = Struct(name='output_data',
                                  mode='cell', data=strain,
                                  dofs=None)
    # Cauchy stress averaged in elements.
    stress = problem.evaluate('ev_cauchy_stress.i1.Omega( solid.D, u )',
                              mode='el_avg')
    out['cauchy_stress'] = Struct(name='output_data',
                                  mode='cell', data=stress,
                                  dofs=None)
    
    return out

# This function will be called by probe.py.
def gen_lines(problem):
    """
    Define three line probes in axial directions.
    
    Parameters
    ----------
    problem : ProblemDefinition instance
        The current ProblemDefinition instance.

    Returns
    -------
    probes : list
        The list of the probes.
    labels : list
        The list of probe labels.
    """
    from sfepy.fem.probes import LineProbe

    mesh = problem.domain.mesh
    bbox = mesh.get_bounding_box()
    cx, cy, cz = 0.5 * bbox.sum(axis=0)
    print bbox
    print cx, cy, cz

    # Probe end points.
    ps0 = [[bbox[0,0], cy, cz],
           [cx, bbox[0,1], cz],
           [cx, cy, bbox[0,2]]]
    ps1 = [[bbox[1,0], cy, cz],
           [cx, bbox[1,1], cz],
           [cx, cy, bbox[1,2]]]
    

    # Use adaptive probe with 10 inital points.
    n_point = -10

    labels = ['%s -> %s' % (p0, p1) for p0, p1 in zip(ps0, ps1)]
    probes = []
    for ip in xrange(len(ps0)):
        p0, p1 = ps0[ip], ps1[ip]
        probes.append(LineProbe(p0, p1, n_point, mesh))

    return probes, labels

# This function will be called by probe.py.
def probe_hook(data, probe, label, problem):
    """
    Parameters
    ----------
    data : dict
        The output data.
    probe : Probe subclass instance
        The probe to be used on data.
    label : str
        The label describing the probe.
    problem : ProblemDefinition instance
        The current ProblemDefinition instance.

    Returns
    -------
    fig : figure
        The matplotlib figure with the probe plot.
    results : dict
        The dict of tuples (pars, vals) of the probe parametrization and the
        corresponding probed data.
    """
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    from sfepy.fem import FieldVariable

    def get_it(name, var_name):
        var = problem.create_variables([var_name])[var_name]
        var.data_from_any(data[name].data)

        pars, vals = probe(var)
        vals = vals.squeeze()
        return pars, vals

    results = {}
    results['u'] = get_it('u', 'u')
    results['cauchy_strain'] = get_it('cauchy_strain', 's')
    results['cauchy_stress'] = get_it('cauchy_stress', 's')

    fig = plt.figure()
    plt.clf()
    fig.subplots_adjust(hspace=0.4)

    plt.subplot(311)
    pars, vals = results['u']
    for ic in range(vals.shape[1]):
        plt.plot(pars, vals[:,ic], label=r'$u_{%d}$' % (ic + 1),
                 lw=1, ls='-', marker='+', ms=3)
    plt.ylabel('displacements')
    plt.xlabel('probe %s' % label, fontsize=8)
    plt.legend(loc='best', prop=fm.FontProperties(size=10))

    sym_indices = ['11', '22', '33', '12', '13', '23']

    plt.subplot(312)
    pars, vals = results['cauchy_strain']
    for ic in range(vals.shape[1]):
        plt.plot(pars, vals[:,ic], label=r'$e_{%s}$' % sym_indices[ic],
                 lw=1, ls='-', marker='+', ms=3)
    plt.ylabel('Cauchy strain')
    plt.xlabel('probe %s' % label, fontsize=8)
    plt.legend(loc='best', prop=fm.FontProperties(size=8))

    plt.subplot(313)
    pars, vals = results['cauchy_stress']
    for ic in range(vals.shape[1]):
        plt.plot(pars, vals[:,ic], label=r'$\tau_{%s}$' % sym_indices[ic],
                 lw=1, ls='-', marker='+', ms=3)
    plt.ylabel('Cauchy stress')
    plt.xlabel('probe %s' % label, fontsize=8)
    plt.legend(loc='best', prop=fm.FontProperties(size=8))

    return plt.gcf(), results

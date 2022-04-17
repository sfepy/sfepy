import pytest
from sfepy.base.testing import check_conditions, run_declaratice_example

examples = [
    'examples/acoustics/acoustics.py',
    'examples/acoustics/acoustics3d.py',
    'examples/acoustics/vibro_acoustic3d.py',
    'examples/diffusion/cube.py',
    'examples/diffusion/darcy_flow_multicomp.py',
    'examples/diffusion/laplace_1d.py',
    'examples/diffusion/laplace_coupling_lcbcs.py',
    'examples/diffusion/laplace_time_ebcs.py',
    'examples/diffusion/poisson.py',
    'examples/diffusion/poisson_field_dependent_material.py',
    'examples/diffusion/poisson_functions.py',
    'examples/diffusion/poisson_neumann.py',
    'examples/diffusion/poisson_periodic_boundary_condition.py',
    'examples/diffusion/sinbc.py',
    'examples/diffusion/time_advection_diffusion.py',
    'examples/diffusion/time_poisson.py',
    'examples/homogenization/linear_elastic_mM.py',
    'examples/large_deformation/active_fibres.py',
    'examples/large_deformation/balloon.py',
    'examples/large_deformation/perfusion_tl.py',
    'examples/linear_elasticity/elastic_contact_planes.py',
    'examples/linear_elasticity/elastic_contact_sphere.py',
    'examples/linear_elasticity/elastic_shifted_periodic.py',
    'examples/linear_elasticity/elastodynamic.py',
    'examples/linear_elasticity/its2D_2.py',
    'examples/linear_elasticity/linear_elastic.py',
    'examples/linear_elasticity/linear_elastic_damping.py',
    'examples/linear_elasticity/linear_elastic_probes.py',
    'examples/linear_elasticity/linear_elastic_tractions.py',
    'examples/linear_elasticity/linear_elastic_up.py',
    'examples/linear_elasticity/linear_viscoelastic.py',
    'examples/linear_elasticity/material_nonlinearity.py',
    'examples/linear_elasticity/nodal_lcbcs.py',
    'examples/linear_elasticity/prestress_fibres.py',
    'examples/linear_elasticity/shell10x_cantilever.py',
    'examples/linear_elasticity/two_bodies_contact.py',
    'examples/multi_physics/biot.py',
    'examples/multi_physics/biot_npbc.py',
    'examples/multi_physics/biot_npbc_lagrange.py',
    'examples/multi_physics/biot_short_syntax.py',
    'examples/multi_physics/piezo_elasticity.py',
    'examples/multi_physics/thermo_elasticity.py',
    'examples/multi_physics/thermo_elasticity_ess.py',
    'examples/navier_stokes/navier_stokes.py',
    'examples/navier_stokes/navier_stokes2d.py',
    'examples/navier_stokes/navier_stokes2d_iga.py',
    'examples/navier_stokes/stabilized_navier_stokes.py',
    'examples/navier_stokes/stokes.py',
    'examples/navier_stokes/stokes_slip_bc.py',
    'examples/quantum/boron.py',
    'examples/quantum/hydrogen.py',
    'examples/quantum/oscillator.py',
    'examples/quantum/well.py',
]

try:
    from igakit import igalib; igalib

except ImportError:
    pass

else:
    examples.extend([
        'examples/diffusion/poisson_iga.py',
        'examples/linear_elasticity/linear_elastic_iga.py',
    ])

@pytest.mark.parametrize('ex_filename', examples)
def test_examples(ex_filename, output_dir):
    conditions = run_declaratice_example(
        ex_filename=ex_filename,
        output_dir=output_dir,
    )
    if ex_filename == 'examples/large_deformation/active_fibres.py':
        # Special-case the first iteration, as the solver converges slowly.
        ok = (conditions[1:] == 0).all()
        ok = ok and (conditions[0] == 1)

    else:
        ok = check_conditions(conditions)

    assert ok

examples_dg = [
    'examples/dg/advection_1D.py',
    'examples/dg/advection_2D.py',
    'examples/dg/advection_diffusion_2D.py',
    'examples/dg/laplace_2D.py',
]
@pytest.mark.parametrize('ex_filename', examples_dg)
def test_examples_dg(ex_filename, output_dir):
    conditions = run_declaratice_example(
        ex_filename=ex_filename,
        output_dir=output_dir,
        ext='.msh',
    )
    ok = check_conditions(conditions)
    assert ok

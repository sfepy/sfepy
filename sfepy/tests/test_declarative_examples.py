import os
import pytest
from sfepy.base.testing import check_conditions, run_declaratice_example

# Relative to sfepy.base_dir, see run_declaratice_example().
examples_dir = 'examples/'

def inedir(filename):
    return os.path.join(examples_dir, filename)

examples = [
    'acoustics/acoustics.py',
    'acoustics/acoustics3d.py',
    'acoustics/helmholtz_apartment.py',
    'acoustics/vibro_acoustic3d.py',
    'diffusion/cube.py',
    'diffusion/darcy_flow_multicomp.py',
    'diffusion/laplace_1d.py',
    'diffusion/laplace_coupling_lcbcs.py',
    'diffusion/laplace_time_ebcs.py',
    'diffusion/poisson.py',
    'diffusion/poisson_field_dependent_material.py',
    'diffusion/poisson_functions.py',
    'diffusion/poisson_neumann.py',
    'diffusion/poisson_periodic_boundary_condition.py',
    'diffusion/sinbc.py',
    'diffusion/time_advection_diffusion.py',
    'diffusion/time_heat_equation_multi_material.py',
    'diffusion/time_poisson.py',
    'homogenization/linear_elastic_mM.py',
    'large_deformation/active_fibres.py',
    'large_deformation/balloon.py',
    'large_deformation/perfusion_tl.py',
    'linear_elasticity/elastic_contact_planes.py',
    'linear_elasticity/elastic_contact_sphere.py',
    'linear_elasticity/elastic_shifted_periodic.py',
    'linear_elasticity/elastodynamic.py',
    'linear_elasticity/its2D_2.py',
    'linear_elasticity/linear_elastic.py',
    'linear_elasticity/linear_elastic_damping.py',
    'linear_elasticity/linear_elastic_probes.py',
    'linear_elasticity/linear_elastic_tractions.py',
    'linear_elasticity/linear_elastic_up.py',
    'linear_elasticity/linear_viscoelastic.py',
    'linear_elasticity/material_nonlinearity.py',
    'linear_elasticity/nodal_lcbcs.py',
    'linear_elasticity/prestress_fibres.py',
    'linear_elasticity/seismic_load.py',
    'linear_elasticity/shell10x_cantilever.py',
    'linear_elasticity/two_bodies_contact.py',
    'multi_physics/biot.py',
    'multi_physics/biot_npbc.py',
    'multi_physics/biot_npbc_lagrange.py',
    'multi_physics/biot_short_syntax.py',
    'multi_physics/piezo_elasticity.py',
    'multi_physics/piezo_elastodynamic.py',
    'multi_physics/thermo_elasticity.py',
    'multi_physics/thermo_elasticity_ess.py',
    'navier_stokes/navier_stokes.py',
    'navier_stokes/navier_stokes2d.py',
    'navier_stokes/stabilized_navier_stokes.py',
    'navier_stokes/stokes.py',
    'navier_stokes/stokes_slip_bc.py',
    'quantum/boron.py',
    'quantum/hydrogen.py',
    'quantum/oscillator.py',
    'quantum/well.py',
]

try:
    from igakit import igalib; igalib

except ImportError:
    pass

else:
    examples.extend([
        'diffusion/poisson_iga.py',
        'linear_elasticity/linear_elastic_iga.py',
        'navier_stokes/navier_stokes2d_iga.py',
    ])

try:
    import primme

except ImportError:
    pass

else:
    examples.extend([
        'linear_elasticity/modal_analysis_declarative.py',
    ])


@pytest.mark.parametrize('ex_filename', examples)
def test_examples(ex_filename, output_dir):
    conditions = run_declaratice_example(
        ex_filename=inedir(ex_filename),
        output_dir=output_dir,
        remove_prefix=examples_dir,
    )
    if ex_filename == 'large_deformation/active_fibres.py':
        # Special-case the first iteration, as the solver converges slowly.
        ok = (conditions[1:] == 0).all()
        ok = ok and (conditions[0] == 1)

    else:
        ok = check_conditions(conditions)

    assert ok

examples_dg = [
    'dg/advection_1D.py',
    'dg/advection_2D.py',
    'dg/advection_diffusion_2D.py',
    'dg/laplace_2D.py',
]
@pytest.mark.parametrize('ex_filename', examples_dg)
def test_examples_dg(ex_filename, output_dir):
    conditions = run_declaratice_example(
        ex_filename=inedir(ex_filename),
        output_dir=output_dir,
        ext='.msh',
        remove_prefix=examples_dir,
    )
    ok = check_conditions(conditions)
    assert ok

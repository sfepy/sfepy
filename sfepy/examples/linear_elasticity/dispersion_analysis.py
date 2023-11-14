#!/usr/bin/env python
"""
Dispersion analysis of a heterogeneous finite scale periodic cell.

The periodic cell mesh has to contain two subdomains Y1 (with the cell ids 1),
Y2 (with the cell ids 2), so that different material properties can be defined
in each of the subdomains (see ``--pars`` option). The command line parameters
can be given in any consistent unit set, for example the basic SI units. The
``--unit-multipliers`` option can be used to rescale the input units to ones
more suitable to the simulation, for example to prevent having different
matrix blocks with large differences of matrix entries magnitudes. The results
are then in the rescaled units.

Usage Examples
--------------

Default material parameters, a square periodic cell with a spherical inclusion,
logs also standard pressure dilatation and shear waves, no eigenvectors::

  python sfepy/examples/linear_elasticity/dispersion_analysis.py meshes/2d/special/circle_in_square.mesh --log-std-waves --eigs-only

As above, with custom eigenvalue solver parameters, and different number of
eigenvalues, mesh size and units used in the calculation::

  python sfepy/examples/linear_elasticity/dispersion_analysis.py meshes/2d/special/circle_in_square.mesh --solver-conf="kind='eig.scipy', method='eigsh', tol=1e-10, maxiter=1000, which='LM', sigma=0" --log-std-waves -n 5 --range=0,640,101 --mode=omega --unit-multipliers=1e-6,1e-2,1e-3 --mesh-size=1e-2 --eigs-only

Default material parameters, a square periodic cell with a square inclusion,
and a very small mesh to allow comparing the omega and kappa modes (full matrix
solver required!)::

  python sfepy/examples/linear_elasticity/dispersion_analysis.py meshes/2d/square_2m.mesh --solver-conf="kind='eig.scipy', method='eigh'" --log-std-waves -n 10 --range=0,640,101 --mesh-size=1e-2 --mode=omega --eigs-only --no-legends --unit-multipliers=1e-6,1e-2,1e-3 -o output/omega

  python sfepy/examples/linear_elasticity/dispersion_analysis.py meshes/2d/square_2m.mesh --solver-conf="kind='eig.qevp', method='companion', mode='inverted', solver={kind='eig.scipy', method='eig'}" --log-std-waves -n 500 --range=0,4000000,1001 --mesh-size=1e-2 --mode=kappa --eigs-only --no-legends --unit-multipliers=1e-6,1e-2,1e-3 -o output/kappa

View/compare the resulting logs::

  python script/plot_logs.py output/omega/frequencies.txt --no-legends -g 1 -o mode-omega.png
  python script/plot_logs.py output/kappa/wave-numbers.txt --no-legends -o mode-kappa.png
  python script/plot_logs.py output/kappa/wave-numbers.txt --no-legends --swap-axes -o mode-kappa-t.png

In contrast to the heterogeneous square periodic cell, a homogeneous
square periodic cell (the region Y2 is empty)::

  python sfepy/examples/linear_elasticity/dispersion_analysis.py meshes/2d/square_1m.mesh --solver-conf="kind='eig.scipy', method='eigh'" --log-std-waves -n 10 --range=0,640,101 --mesh-size=1e-2 --mode=omega --eigs-only --no-legends --unit-multipliers=1e-6,1e-2,1e-3 -o output/omega-h

  python script/plot_logs.py output/omega-h/frequencies.txt --no-legends -g 1 -o mode-omega-h.png

Use the Brillouin stepper::

  python sfepy/examples/linear_elasticity/dispersion_analysis.py meshes/2d/special/circle_in_square.mesh --log-std-waves -n=60 --eigs-only --no-legends --stepper=brillouin

  python script/plot_logs.py output/frequencies.txt -g 0 --rc="'font.size':14, 'lines.linewidth' : 3, 'lines.markersize' : 4" -o brillouin-stepper-kappas.png

  python script/plot_logs.py output/frequencies.txt -g 1 --no-legends --rc="'font.size':14, 'lines.linewidth' : 3, 'lines.markersize' : 4" -o brillouin-stepper-omegas.png

Additional arguments can be passed to the problem configuration's
:func:`define()` function using the ``--define-kwargs`` option. In this file,
only the mesh vertex separation parameter `mesh_eps` can be used::

  python sfepy/examples/linear_elasticity/dispersion_analysis.py meshes/2d/special/circle_in_square.mesh --log-std-waves --eigs-only --define-kwargs="mesh_eps=1e-10" --save-regions
"""
from __future__ import absolute_import
import os
import sys
sys.path.append('.')
import gc
from copy import copy
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import numpy as nm
import matplotlib.pyplot as plt

from sfepy.base.base import import_file, output, Struct
from sfepy.base.conf import dict_from_string, ProblemConf
from sfepy.base.ioutils import ensure_path, remove_files_patterns, save_options
from sfepy.base.log import Log
from sfepy.discrete.fem import MeshIO
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson as stiffness
import sfepy.mechanics.matcoefs as mc
from sfepy.mechanics.units import apply_unit_multipliers, apply_units_to_pars
import sfepy.discrete.fem.periodic as per
from sfepy.discrete.fem.meshio import convert_complex_output
from sfepy.homogenization.utils import define_box_regions
from sfepy.discrete import Problem
from sfepy.mechanics.tensors import get_von_mises_stress
from sfepy.solvers import Solver
from sfepy.solvers.ts import get_print_info, TimeStepper
from sfepy.linalg.utils import output_array_stats, max_diff_csr

pars_kinds = {
    'young1' : 'stress',
    'poisson1' : 'one',
    'density1' : 'density',
    'young2' : 'stress',
    'poisson2' : 'one',
    'density2' : 'density',
}

def define(filename_mesh, pars, approx_order, refinement_level, solver_conf,
           plane='strain', post_process=False, mesh_eps=1e-8):
    io = MeshIO.any_from_filename(filename_mesh)
    bbox = io.read_bounding_box()
    dim = bbox.shape[1]

    options = {
        'absolute_mesh_path' : True,
        'refinement_level' : refinement_level,
        'allow_empty_regions' : True,
        'post_process_hook' : 'compute_von_mises' if post_process else None,
    }

    fields = {
        'displacement': ('complex', dim, 'Omega', approx_order),
    }

    materials = {
        'm' : ({
            'D' : {'Y1' : stiffness(dim,
                                    young=pars.young1,
                                    poisson=pars.poisson1,
                                    plane=plane),
                   'Y2' : stiffness(dim,
                                    young=pars.young2,
                                    poisson=pars.poisson2,
                                    plane=plane)},
            'density' : {'Y1' : pars.density1, 'Y2' : pars.density2},
        },),
        'wave' : 'get_wdir',
    }

    variables = {
        'u' : ('unknown field', 'displacement', 0),
        'v' : ('test field', 'displacement', 'u'),
    }

    regions = {
        'Omega' : 'all',
        'Y1': 'cells of group 1',
        'Y2': 'cells of group 2',
    }
    regions.update(define_box_regions(dim,
                                      bbox[0], bbox[1], mesh_eps))

    ebcs = {
    }

    if dim == 3:
        epbcs = {
            'periodic_x' : (['Left', 'Right'], {'u.all' : 'u.all'},
                            'match_x_plane'),
            'periodic_y' : (['Near', 'Far'], {'u.all' : 'u.all'},
                            'match_y_plane'),
            'periodic_z' : (['Top', 'Bottom'], {'u.all' : 'u.all'},
                            'match_z_plane'),
        }
    else:
        epbcs = {
            'periodic_x' : (['Left', 'Right'], {'u.all' : 'u.all'},
                            'match_y_line'),
            'periodic_y' : (['Bottom', 'Top'], {'u.all' : 'u.all'},
                            'match_x_line'),
        }

    per.set_accuracy(mesh_eps)
    functions = {
        'match_x_plane' : (per.match_x_plane,),
        'match_y_plane' : (per.match_y_plane,),
        'match_z_plane' : (per.match_z_plane,),
        'match_x_line' : (per.match_x_line,),
        'match_y_line' : (per.match_y_line,),
        'get_wdir' : (get_wdir,),
    }

    integrals = {
        'i' : 2 * approx_order,
    }

    equations = {
        'K' : 'dw_lin_elastic.i.Omega(m.D, v, u)',
        'S' : 'dw_elastic_wave.i.Omega(m.D, wave.vec, v, u)',
        'R' : """1j * dw_elastic_wave_cauchy.i.Omega(m.D, wave.vec, u, v)
               - 1j * dw_elastic_wave_cauchy.i.Omega(m.D, wave.vec, v, u)""",
        'M' : 'dw_dot.i.Omega(m.density, v, u)',
    }

    solver_0 = solver_conf.copy()
    solver_0['name'] = 'eig'

    return locals()

def get_wdir(ts, coors, mode=None,
             equations=None, term=None, problem=None, wdir=None, **kwargs):
    if mode == 'special':
        return {'vec' : wdir}

def set_wave_dir(pb, wdir):
    materials = pb.get_materials()
    wave_mat = materials['wave']
    wave_mat.set_extra_args(wdir=wdir)

def save_materials(output_dir, pb, options):
    stiffness = pb.evaluate('ev_integrate_mat.2.Omega(m.D, u)',
                            mode='el_avg', copy_materials=False, verbose=False)
    young, poisson = mc.youngpoisson_from_stiffness(stiffness,
                                                    plane=options.plane)
    density = pb.evaluate('ev_integrate_mat.2.Omega(m.density, u)',
                          mode='el_avg', copy_materials=False, verbose=False)

    out = {}
    out['young'] = Struct(name='young', mode='cell',
                          data=young[..., None, None])
    out['poisson'] = Struct(name='poisson', mode='cell',
                            data=poisson[..., None, None])
    out['density'] = Struct(name='density', mode='cell', data=density)
    materials_filename = os.path.join(output_dir, 'materials.vtk')
    pb.save_state(materials_filename, out=out)

def get_std_wave_fun(pb, options):
    stiffness = pb.evaluate('ev_integrate_mat.2.Omega(m.D, u)',
                            mode='el_avg', copy_materials=False, verbose=False)
    young, poisson = mc.youngpoisson_from_stiffness(stiffness,
                                                    plane=options.plane)
    density = pb.evaluate('ev_integrate_mat.2.Omega(m.density, u)',
                          mode='el_avg', copy_materials=False, verbose=False)

    lam, mu = mc.lame_from_youngpoisson(young, poisson,
                                        plane=options.plane)
    alam = nm.average(lam)
    amu = nm.average(mu)
    adensity = nm.average(density)

    cp = nm.sqrt((alam + 2.0 * amu) / adensity)
    cs = nm.sqrt(amu / adensity)
    output('average p-wave speed:', cp)
    output('average shear wave speed:', cs)

    log_names = [r'$\omega_p$', r'$\omega_s$']
    log_plot_kwargs = [{'ls' : '--', 'color' : 'k'},
                       {'ls' : '--', 'color' : 'gray'}]

    if options.mode == 'omega':
        fun = lambda wmag, wdir: (cp * wmag, cs * wmag)

    else:
        fun = lambda wmag, wdir: (wmag / cp, wmag / cs)

    return fun, log_names, log_plot_kwargs

def get_stepper(rng, pb, options):
    if options.stepper == 'linear':
        stepper = TimeStepper(rng[0], rng[1], dt=None, n_step=rng[2])
        return stepper

    bbox = pb.domain.mesh.get_bounding_box()

    bzone = 2.0 * nm.pi / (bbox[1] - bbox[0])

    num = rng[2] // 3

    class BrillouinStepper(Struct):
        """
        Step over 1. Brillouin zone in xy plane.
        """
        def __init__(self, t0, t1, dt=None, n_step=None, step=None, **kwargs):
            Struct.__init__(self, t0=t0, t1=t1, dt=dt, n_step=n_step, step=step)

            self.n_digit, self.format, self.suffix = get_print_info(self.n_step)

        def __iter__(self):
            ts = TimeStepper(0, bzone[0], dt=None, n_step=num)
            for ii, val in ts:
                yield ii, val, nm.array([1.0, 0.0])
                if ii == (num-2): break

            ts = TimeStepper(0, bzone[1], dt=None, n_step=num)
            for ii, k1 in ts:
                wdir = nm.array([bzone[0], k1])
                val = nm.linalg.norm(wdir)
                wdir = wdir / val
                yield num + ii, val, wdir
                if ii == (num-2): break

            wdir = nm.array([bzone[0], bzone[1]])
            val = nm.linalg.norm(wdir)
            wdir = wdir / val
            ts = TimeStepper(0, 1, dt=None, n_step=num)
            for ii, _ in ts:
                yield 2 * num + ii, val * (1.0 - float(ii)/(num-1)), wdir

    stepper = BrillouinStepper(0, 1, n_step=rng[2])

    return stepper

def compute_von_mises(out, pb, state, extend=False, wmag=None, wdir=None):
    """
    Calculate the von Mises stress.
    """
    stress = pb.evaluate('ev_cauchy_stress.i.Omega(m.D, u)', mode='el_avg')

    vms = get_von_mises_stress(stress.squeeze())
    vms.shape = (vms.shape[0], 1, 1, 1)
    out['von_mises_stress'] = Struct(name='output_data', mode='cell',
                                     data=vms)

    return out

def save_eigenvectors(filename, svecs, wmag, wdir, pb):
    if svecs is None: return

    variables = pb.set_default_state()
    # Make full eigenvectors (add DOFs fixed by boundary conditions).
    vecs = nm.empty((variables.di.n_dof_total, svecs.shape[1]),
                    dtype=svecs.dtype)
    for ii in range(svecs.shape[1]):
        vecs[:, ii] = variables.make_full_vec(svecs[:, ii])

    # Save the eigenvectors.
    out = {}

    pp_name = pb.conf.options.get('post_process_hook')
    pp = getattr(pb.conf.funmod, pp_name if pp_name is not None else '',
                 lambda out, *args, **kwargs: out)

    for ii in range(svecs.shape[1]):
        variables.set_state(vecs[:, ii])
        aux = variables.create_output()
        aux2 = {}
        pp(aux2, pb, variables, wmag=wmag, wdir=wdir)
        aux.update(convert_complex_output(aux2))
        out.update({key + '%03d' % ii : aux[key] for key in aux})

    pb.save_state(filename, out=out)

def assemble_matrices(define, mod, pars, set_wave_dir, options, wdir=None):
    """
    Assemble the blocks of dispersion eigenvalue problem matrices.
    """
    define_dict = define(filename_mesh=options.mesh_filename,
                         pars=pars,
                         approx_order=options.order,
                         refinement_level=options.refine,
                         solver_conf=options.solver_conf,
                         plane=options.plane,
                         post_process=options.post_process,
                         **options.define_kwargs)

    conf = ProblemConf.from_dict(define_dict, mod)

    pb = Problem.from_conf(conf)
    pb.dispersion_options = options
    pb.set_output_dir(options.output_dir)
    dim = pb.domain.shape.dim

    # Set the normalized wave vector direction to the material(s).
    if wdir is None:
        wdir = nm.asarray(options.wave_dir[:dim], dtype=nm.float64)
        wdir = wdir / nm.linalg.norm(wdir)
    set_wave_dir(pb, wdir)

    bbox = pb.domain.mesh.get_bounding_box()
    size = (bbox[1] - bbox[0]).max()
    scaling0 = apply_unit_multipliers([1.0], ['length'],
                                      options.unit_multipliers)[0]
    scaling = scaling0
    if options.mesh_size is not None:
        scaling *= options.mesh_size / size
    output('scaling factor of periodic cell mesh coordinates:', scaling)
    output('new mesh size with applied unit multipliers:', scaling * size)
    pb.domain.mesh.coors[:] *= scaling
    pb.set_mesh_coors(pb.domain.mesh.coors, update_fields=True)

    bzone = 2.0 * nm.pi / (scaling * size)
    output('1. Brillouin zone size:', bzone * scaling0)
    output('1. Brillouin zone size with applied unit multipliers:', bzone)

    pb.time_update()
    pb.update_materials()

    # Assemble the matrices.
    mtxs = {}
    for key, eq in pb.equations.iteritems():
        mtxs[key] = mtx = pb.mtx_a.copy()
        mtx = eq.evaluate(mode='weak', dw_mode='matrix', asm_obj=mtx)
        mtx.eliminate_zeros()
        output_array_stats(mtx.data, 'nonzeros in %s' % key)

        output('symmetry checks:')
        output('%s - %s^T:' % (key, key), max_diff_csr(mtx, mtx.T))
        output('%s - %s^H:' % (key, key), max_diff_csr(mtx, mtx.H))

    return pb, wdir, bzone, mtxs

def setup_n_eigs(options, pb, mtxs):
    """
    Setup the numbers of eigenvalues based on options and numbers of DOFs.
    """
    solver_n_eigs = n_eigs = options.n_eigs
    n_dof = mtxs['K'].shape[0]
    if options.mode == 'omega':
        if options.n_eigs > n_dof:
            n_eigs = n_dof
            solver_n_eigs = None

    else:
        if options.n_eigs > 2 * n_dof:
            n_eigs = 2 * n_dof
            solver_n_eigs = None

    return solver_n_eigs, n_eigs

def build_evp_matrices(mtxs, val, mode, pb):
    """
    Build the matrices of the dispersion eigenvalue problem.
    """
    if mode == 'omega':
        mtx_a = mtxs['K'] + val**2 * mtxs['S'] + val * mtxs['R']
        output('A - A^H:', max_diff_csr(mtx_a, mtx_a.H))

        evp_mtxs = (mtx_a, mtxs['M'])

    else:
        evp_mtxs = (mtxs['S'], mtxs['R'], mtxs['K'] - val**2 * mtxs['M'])

    return evp_mtxs

def process_evp_results(eigs, svecs, val, wdir, bzone, pb, mtxs, options,
                        std_wave_fun=None):
    """
    Transform eigenvalues to either omegas or kappas, depending on `mode`.
    Transform eigenvectors, if available, depending on `mode`.
    Return also the values to log.
    """
    if options.mode == 'omega':
        omegas = nm.sqrt(eigs)

        output('eigs, omegas:')
        for ii, om in enumerate(omegas):
            output('{:>3}. {: .10e}, {:.10e}'.format(ii, eigs[ii], om))

        if options.stepper == 'linear':
            out = tuple(eigs) + tuple(omegas)

        else:
            out = tuple(val * wdir) + tuple(omegas)

        if std_wave_fun is not None:
            out = out + std_wave_fun(val, wdir)

        return omegas, svecs, out

    else:
        kappas = eigs.copy()
        rks = kappas.copy()

        # Mask modes far from 1. Brillouin zone.
        max_kappa = 1.2 * bzone
        kappas[kappas.real > max_kappa] = nm.nan

        # Mask non-physical modes.
        kappas[kappas.real < 0] = nm.nan
        kappas[nm.abs(kappas.imag) > 1e-10] = nm.nan
        out = tuple(kappas.real)

        output('raw kappas, masked real part:',)
        for ii, kr in enumerate(kappas.real):
            output('{:>3}. {: 23.5e}, {:.10e}'.format(ii, rks[ii], kr))

        if svecs is not None:
            n_dof = mtxs['K'].shape[0]
            # Select only vectors corresponding to physical modes.
            ii = nm.isfinite(kappas.real)
            svecs = svecs[:n_dof, ii]

        if std_wave_fun is not None:
            out = out + tuple(ii if ii <= max_kappa else nm.nan
                              for ii in std_wave_fun(val, wdir))

        return kappas, svecs, out

helps = {
    'pars' :
    'material parameters in Y1, Y2 subdomains in basic units.'
    ' The default parameters are:'
    ' young1, poisson1, density1, young2, poisson2, density2'
    ' [default: %(default)s]',
    'conf' :
    'if given, an alternative problem description file with apply_units() and'
    ' define() functions [default: %(default)s]',
    'define_kwargs' : 'additional keyword arguments passed to define()',
    'mesh_size' :
    'desired mesh size (max. of bounding box dimensions) in basic units'
    ' - the input periodic cell mesh is rescaled to this size'
    ' [default: %(default)s]',
    'unit_multipliers' :
    'basic unit multipliers (time, length, mass) [default: %(default)s]',
    'plane' :
    'for 2D problems, plane strain or stress hypothesis selection'
    ' [default: %(default)s]',
    'wave_dir' : 'the wave vector direction (will be normalized)'
    ' [default: %(default)s]',
    'mode' : 'solution mode: omega = solve a generalized EVP for omega,'
    ' kappa = solve a quadratic generalized EVP for kappa'
    ' [default: %(default)s]',
    'stepper' : 'the range stepper. For "brillouin", only the number'
    ' of items from --range is used'
    ' [default: %(default)s]',
    'range' : 'the wave vector magnitude / frequency range'
    ' (like numpy.linspace) depending on the mode option'
    ' [default: %(default)s]',
    'order' : 'displacement field approximation order [default: %(default)s]',
    'refine' : 'number of uniform mesh refinements [default: %(default)s]',
    'n_eigs' : 'the number of eigenvalues to compute [default: %(default)s]',
    'eigs_only' : 'compute only eigenvalues, not eigenvectors',
    'post_process' : 'post-process eigenvectors',
    'solver_conf' : 'eigenvalue problem solver configuration options'
    ' [default: %(default)s]',
    'save_regions' : 'save defined regions into'
    ' <output_directory>/regions.vtk',
    'save_materials' : 'save material parameters into'
    ' <output_directory>/materials.vtk',
    'log_std_waves' : 'log also standard pressure dilatation and shear waves',
    'no_legends' :
    'do not show legends in the log plots',
    'no_show' :
    'do not show the log figure',
    'silent' : 'do not print messages to screen',
    'clear' :
    'clear old solution files from output directory',
    'output_dir' :
    'output directory [default: %(default)s]',
    'mesh_filename' :
    'input periodic cell mesh file name [default: %(default)s]',
}

def main():
    # Aluminium and epoxy.
    default_pars = '70e9,0.35,2.799e3,3.8e9,0.27,1.142e3'
    default_solver_conf = ("kind='eig.scipy',method='eigsh',tol=1.0e-5,"
                           "maxiter=1000,which='LM',sigma=0.0")

    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--pars', metavar='name1=value1,name2=value2,...'
                        ' or value1,value2,...',
                        action='store', dest='pars',
                        default=default_pars, help=helps['pars'])
    parser.add_argument('--conf', metavar='filename',
                        action='store', dest='conf',
                        default=None, help=helps['conf'])
    parser.add_argument('--define-kwargs', metavar='dict-like',
                        action='store', dest='define_kwargs',
                        default=None, help=helps['define_kwargs'])
    parser.add_argument('--mesh-size', type=float, metavar='float',
                        action='store', dest='mesh_size',
                        default=None, help=helps['mesh_size'])
    parser.add_argument('--unit-multipliers',
                        metavar='c_time,c_length,c_mass',
                        action='store', dest='unit_multipliers',
                        default='1.0,1.0,1.0', help=helps['unit_multipliers'])
    parser.add_argument('--plane', action='store', dest='plane',
                        choices=['strain', 'stress'],
                        default='strain', help=helps['plane'])
    parser.add_argument('--wave-dir', metavar='float,float[,float]',
                        action='store', dest='wave_dir',
                        default='1.0,0.0,0.0', help=helps['wave_dir'])
    parser.add_argument('--mode', action='store', dest='mode',
                        choices=['omega', 'kappa'],
                        default='omega', help=helps['mode'])
    parser.add_argument('--stepper', action='store', dest='stepper',
                        choices=['linear', 'brillouin'],
                        default='linear', help=helps['stepper'])
    parser.add_argument('--range', metavar='start,stop,count',
                        action='store', dest='range',
                        default='0,6.4,33', help=helps['range'])
    parser.add_argument('--order', metavar='int', type=int,
                        action='store', dest='order',
                        default=1, help=helps['order'])
    parser.add_argument('--refine', metavar='int', type=int,
                        action='store', dest='refine',
                        default=0, help=helps['refine'])
    parser.add_argument('-n', '--n-eigs', metavar='int', type=int,
                        action='store', dest='n_eigs',
                        default=6, help=helps['n_eigs'])
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--eigs-only',
                       action='store_true', dest='eigs_only',
                       default=False, help=helps['eigs_only'])
    group.add_argument('--post-process',
                       action='store_true', dest='post_process',
                       default=False, help=helps['post_process'])
    parser.add_argument('--solver-conf', metavar='dict-like',
                        action='store', dest='solver_conf',
                        default=default_solver_conf, help=helps['solver_conf'])
    parser.add_argument('--save-regions',
                        action='store_true', dest='save_regions',
                        default=False, help=helps['save_regions'])
    parser.add_argument('--save-materials',
                        action='store_true', dest='save_materials',
                        default=False, help=helps['save_materials'])
    parser.add_argument('--log-std-waves',
                        action='store_true', dest='log_std_waves',
                        default=False, help=helps['log_std_waves'])
    parser.add_argument('--no-legends',
                        action='store_false', dest='show_legends',
                        default=True, help=helps['no_legends'])
    parser.add_argument('--no-show',
                        action='store_false', dest='show',
                        default=True, help=helps['no_show'])
    parser.add_argument('--silent',
                        action='store_true', dest='silent',
                        default=False, help=helps['silent'])
    parser.add_argument('-c', '--clear',
                        action='store_true', dest='clear',
                        default=False, help=helps['clear'])
    parser.add_argument('-o', '--output-dir', metavar='path',
                        action='store', dest='output_dir',
                        default='output', help=helps['output_dir'])
    parser.add_argument('mesh_filename', default='',
                        help=helps['mesh_filename'])
    options = parser.parse_args()

    output_dir = options.output_dir

    output.set_output(filename=os.path.join(output_dir,'output_log.txt'),
                      combined=options.silent == False)

    if options.conf is not None:
        mod = import_file(options.conf)

    else:
        mod = sys.modules[__name__]

    pars_kinds = mod.pars_kinds
    define = mod.define
    set_wave_dir = mod.set_wave_dir
    setup_n_eigs = mod.setup_n_eigs
    build_evp_matrices = mod.build_evp_matrices
    save_materials = mod.save_materials
    get_std_wave_fun = mod.get_std_wave_fun
    get_stepper = mod.get_stepper
    process_evp_results = mod.process_evp_results
    save_eigenvectors = mod.save_eigenvectors

    try:
        options.pars = dict_from_string(options.pars)

    except:
        aux = [float(ii) for ii in options.pars.split(',')]
        options.pars = {key : aux[ii]
                        for ii, key in enumerate(pars_kinds.keys())}

    options.unit_multipliers = [float(ii)
                                for ii in options.unit_multipliers.split(',')]
    options.wave_dir = [float(ii)
                        for ii in options.wave_dir.split(',')]
    aux = options.range.split(',')
    options.range = [float(aux[0]), float(aux[1]), int(aux[2])]
    options.solver_conf = dict_from_string(options.solver_conf)
    options.define_kwargs = dict_from_string(options.define_kwargs)

    if options.clear:
        remove_files_patterns(output_dir,
                              ['*.h5', '*.vtk', '*.txt'],
                              ignores=['output_log.txt'],
                              verbose=True)

    filename = os.path.join(output_dir, 'options.txt')
    ensure_path(filename)
    save_options(filename, [('options', vars(options))],
                 quote_command_line=True)

    pars = apply_units_to_pars(options.pars, pars_kinds,
                               options.unit_multipliers)
    output('material parameter names and kinds:')
    output(pars_kinds)
    output('material parameters with applied unit multipliers:')
    output(pars)

    pars = Struct(**pars)

    if options.mode == 'omega':
        rng = copy(options.range)
        rng[:2] = apply_unit_multipliers(options.range[:2],
                                         ['wave_number', 'wave_number'],
                                         options.unit_multipliers)
        output('wave number range with applied unit multipliers:', rng)

    else:
        if options.stepper == 'brillouin':
            raise ValueError('Cannot use "brillouin" stepper in kappa mode!')

        rng = copy(options.range)
        rng[:2] = apply_unit_multipliers(options.range[:2],
                                         ['frequency', 'frequency'],
                                         options.unit_multipliers)
        output('frequency range with applied unit multipliers:', rng)

    pb, wdir, bzone, mtxs = assemble_matrices(define, mod, pars, set_wave_dir,
                                              options)
    dim = pb.domain.shape.dim

    if dim != 2:
        options.plane = 'strain'

    if options.save_regions:
        pb.save_regions_as_groups(os.path.join(output_dir, 'regions'))

    if options.save_materials:
        save_materials(output_dir, pb, options)

    conf = pb.solver_confs['eig']
    eig_solver = Solver.any_from_conf(conf)

    n_eigs, options.n_eigs = setup_n_eigs(options, pb, mtxs)

    get_color = lambda ii: plt.cm.viridis((float(ii)
                                           / (max(options.n_eigs, 2) - 1)))
    plot_kwargs = [{'color' : get_color(ii), 'ls' : '', 'marker' : 'o'}
                  for ii in range(options.n_eigs)]
    get_color_dim = lambda ii: plt.cm.viridis((float(ii) / (max(dim, 2) -1)))
    plot_kwargs_dim = [{'color' : get_color_dim(ii), 'ls' : '', 'marker' : 'o'}
                       for ii in range(dim)]

    log_names = []
    log_plot_kwargs = []
    if options.log_std_waves:
        std_wave_fun, log_names, log_plot_kwargs = get_std_wave_fun(
            pb, options)

    else:
        std_wave_fun = None

    stepper = get_stepper(rng, pb, options)

    if options.mode == 'omega':
        eigenshapes_filename = os.path.join(output_dir,
                                            'frequency-eigenshapes-%s.vtk'
                                            % stepper.suffix)

        if options.stepper == 'linear':
            log = Log([[r'$\lambda_{%d}$' % ii for ii in range(options.n_eigs)],
                   [r'$\omega_{%d}$'
                    % ii for ii in range(options.n_eigs)] + log_names],
                  plot_kwargs=[plot_kwargs, plot_kwargs + log_plot_kwargs],
                  formats=[['{:.12e}'] * options.n_eigs,
                           ['{:.12e}'] * (options.n_eigs + len(log_names))],
                  yscales=['linear', 'linear'],
                  xlabels=[r'$\kappa$', r'$\kappa$'],
                  ylabels=[r'eigenvalues $\lambda_i$',
                           r'frequencies $\omega_i$'],
                  show_legends=options.show_legends,
                  is_plot=options.show,
                  log_filename=os.path.join(output_dir, 'frequencies.txt'),
                  aggregate=1000, sleep=0.1)

        else:
            log = Log([[r'$\kappa_{%d}$'% ii for ii in range(dim)],
                       [r'$\omega_{%d}$'
                        % ii for ii in range(options.n_eigs)] + log_names],
                      plot_kwargs=[plot_kwargs_dim,
                                   plot_kwargs + log_plot_kwargs],
                      formats=[['{:.12e}'] * dim,
                               ['{:.12e}'] * (options.n_eigs + len(log_names))],
                      yscales=['linear', 'linear'],
                      xlabels=[r'', r''],
                      ylabels=[r'wave vector $\kappa$',
                               r'frequencies $\omega_i$'],
                      show_legends=options.show_legends,
                      is_plot=options.show,
                      log_filename=os.path.join(output_dir, 'frequencies.txt'),
                      aggregate=1000, sleep=0.1)

        for aux in stepper:
            if options.stepper == 'linear':
                iv, wmag = aux

            else:
                iv, wmag, wdir = aux

            output('step %d: wave vector %s' % (iv, wmag * wdir))

            if options.stepper == 'brillouin':
                pb, _, bzone, mtxs = assemble_matrices(
                    define, mod, pars, set_wave_dir, options, wdir=wdir)

            evp_mtxs = build_evp_matrices(mtxs, wmag, options.mode, pb)

            if options.eigs_only:
                eigs = eig_solver(*evp_mtxs, n_eigs=n_eigs,
                                  eigenvectors=False)
                svecs = None

            else:
                eigs, svecs = eig_solver(*evp_mtxs, n_eigs=n_eigs,
                                         eigenvectors=True)

            omegas, svecs, out = process_evp_results(
                eigs, svecs, wmag, wdir, bzone, pb, mtxs, options,
                std_wave_fun=std_wave_fun
            )
            if options.stepper == 'linear':
                log(*out, x=[wmag, wmag])

            else:
                log(*out, x=[iv, iv])

            save_eigenvectors(eigenshapes_filename % iv, svecs, wmag, wdir, pb)

            gc.collect()

        log(save_figure=os.path.join(output_dir, 'frequencies.png'))
        log(finished=True)

    else:
        eigenshapes_filename = os.path.join(output_dir,
                                            'wave-number-eigenshapes-%s.vtk'
                                            % stepper.suffix)

        log = Log([[r'$\kappa_{%d}$' % ii for ii in range(options.n_eigs)]
                   + log_names],
                  plot_kwargs=[plot_kwargs + log_plot_kwargs],
                  formats=[['{:.12e}'] * (options.n_eigs + len(log_names))],
                  yscales=['linear'],
                  xlabels=[r'$\omega$'],
                  ylabels=[r'wave numbers $\kappa_i$'],
                  show_legends=options.show_legends,
                  is_plot=options.show,
                  log_filename=os.path.join(output_dir, 'wave-numbers.txt'),
                  aggregate=1000, sleep=0.1)
        for io, omega in stepper:
            output('step %d: frequency %s' % (io, omega))

            evp_mtxs = build_evp_matrices(mtxs, omega, options.mode, pb)

            if options.eigs_only:
                eigs = eig_solver(*evp_mtxs, n_eigs=n_eigs,
                                  eigenvectors=False)
                svecs = None

            else:
                eigs, svecs = eig_solver(*evp_mtxs, n_eigs=n_eigs,
                                         eigenvectors=True)

            kappas, svecs, out = process_evp_results(
                eigs, svecs, omega, wdir, bzone, pb, mtxs, options,
                std_wave_fun=std_wave_fun
            )
            log(*out, x=[omega])

            save_eigenvectors(eigenshapes_filename % io, svecs, kappas, wdir,
                              pb)

            gc.collect()

        log(save_figure=os.path.join(output_dir, 'wave-numbers.png'))
        log(finished=True)

if __name__ == '__main__':
    main()

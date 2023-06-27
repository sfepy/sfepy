#!/usr/bin/env python
"""
Generate the images and rst files for gallery of SfePy examples.

The following steps need to be made to regenerate the documentation with the
updated example files:

1. remove doc/examples/*::

   $ rm -rf doc/examples/*

2. generate the files:

   $ ./tools/gen_gallery.py

3. regenerate the documentation::

   $ python setup.py htmldocs
"""
from __future__ import absolute_import
import sys
import six
sys.path.append('.')
import os
import tempfile
import glob
import re
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import numpy as nm
import matplotlib.image as image

import sfepy
from sfepy.base.base import (get_default, ordered_iteritems,
                             import_file, output, Struct)
from sfepy.base.ioutils import (ensure_path, locate_files, remove_files,
                                edit_filename)
from sfepy.scripts.resview import pv_plot, get_camera_position
import pyvista as pv

omits = [
    '__init__.py',
]

omit_images = [
    'vibro_acoustic3d_mid.py',
    'its2D_5.py',
    'linear_elastic_probes.py',
]

omit_dirs = [
    re.compile('.*output.*/').match,
]

custom = {
    'acoustics/vibro_acoustic3d.py': {
        '_Gamma0_1': {'view_2d': True, 'max_plots': 2},
        '_Omega1': {'camera': [45, 55, 0.8]},
        '_Omega2': {'camera': [45, 55, 0.8]},
    },
    'acoustics/acoustics3d.py': {
        '_Omega_1': {
            'camera': [75, 135, 1.4],
            'grid_vector1': [1.2, 0, 0],
        },
        '_Omega_2': {
            'camera': [75, 135, 1],
            'grid_vector1': [1.2, 0, 0],
        },
    },
    'acoustics/helmholtz_apartment.py': {
        '': {
            'fields': ['imag.E:wimag.E:f10%:p0', 'mat_id:p1'],
            'force_view_3d': True,
            'grid_vector1': [1.1, 0, 0],
            'camera_position': [-9.22684,-8.37688,10.7623,
                                1.51653,0.122742,-0.646791,
                                0.472758,0.432784,0.767593],
            'color_map': 'seismic',
        },
    },
    'diffusion/cube.py': {
        '': {'camera': [225, 55, 0.7]},
    },
    'diffusion/laplace_time_ebcs.py': {
        '': {'camera': [225, 55, 0.7]},
    },
    'diffusion/laplace_fluid_2d.py': {
        '': {'fields': ['phi:p0', 'phi:t50:p0']},
    },
    'diffusion/laplace_1d.py': {
        '': {'fields': ['t:wt:p0', 't:p0']},
    },
    'diffusion/laplace_coupling_lcbcs.py': {
        '': {
            'fields': ['u1:wu1:p0', 'u1:vw:p0',
                       'u2:wu2:p1', 'u2:vw:p1'],
            'force_view_3d': True,
        },
    },
    'diffusion/poisson_iga.py': {
        '': {
            'fields': ['t:wt:p0', 't:vw:o0.4:p0'],
            'force_view_3d': True,
            'camera': [30, 60, 1.],
        },
    },
    'diffusion/sinbc.py': {
        '_t': {
            'fields': ['t:wt:p0', '1:vw:p0'],
            'force_view_3d': True,
            'camera': [0, -45, 1.],
        },
        '_grad': {
            'fields': ['grad:g:f0.01:p0', '1:vw:p0'],
            'force_view_3d': True,
            'camera': [0, -45, 1.5],
        },
    },
    'diffusion/time_heat_equation_multi_material.py': {
        '': {
            'isosurfaces': 10,
            'outline': True,
            'camera': [-50, -230, 1],
        },
    },
    'linear_elasticity/elastic_contact_planes.py': {
        '': {
            'fields': ['u:wu:p0', '1:vw:p0'],
            'camera': [225, 55, 0.7],
        },
    },
    'linear_elasticity/elastic_contact_sphere.py': {
        '': {
            'fields': ['u:wu:p0', '1:vw:p0'],
            'camera': [225, 55, 0.7],
        },
    },
    'linear_elasticity/elastic_shifted_periodic.py': {
        '': {'fields': ['von_mises_stress:r:wu:p0', '1:vw:p0']},
    },
    'linear_elasticity/elastodynamic.py': {
        '': {
            'fields': ['u:wu:f1e3:p0', '1:vw:p0',
                       'cauchy_strain:p1', 'cauchy_stress:p2'],
        },
    },
    'linear_elasticity/linear_elastic_iga.py': {
        '': {
            'fields': ['u:wu:p0', '1:vw:p0'],
            'camera': [-45, 55, 1],
        },
    },
    'linear_elasticity/linear_viscoelastic.py': {
        '': {'camera': [225, 75, 0.88]}
    },
    'linear_elasticity/modal_analysis_declarative.py': {
        '': {
            'fields': ['u003:wu003:f30%:p0', '1:vw:p0'],
            'camera_position': [-2.30562,-2.2604,-0.325838,
                                -0.0714771,0.0374911,0.0287214,
                                -0.0245554,-0.129086,0.991329],
        },
    },
    'linear_elasticity/seismic_load.py': {
        '': {
            'fields': ['cauchy_stress:wu:f10:p0', '1:vw:p0'],
        },
    },
    'linear_elasticity/shell10x_cantilever.py': {
        '': {
            'fields': ['u_disp:wu_disp:p0', '1:vw:p0',
                       'u_rot:p1', '1:vw:p1'],
            'camera': [-45, 75, 1],
            'grid_vector1': [1, 0, 0],
        },
    },
    'navier_stokes/stokes_slip_bc.py': {
        '': {
            'fields': ['u:g:f.25:p0', 'u:o.4:p0', 'p:p1'],
            'camera': [-45, 55, 1],
            'grid_vector1': [0, 1.2, 0]
        },
    },
    'multi_physics/piezo_elasticity.py': {
        '': {'fields': ['u:g:p0', 'cauchy_strain:p1',
                        'elastic_stress:p2', 'piezo_stress:p3'],
             'grid_vector1': [1.2, 0, 0],
             'max_plots': 4},
    },
    'multi_physics/piezo_elastodynamic.py': {
        '': {
            'fields': ['p:wu:f2000:p0', '1:vw:wu:f2000:p0'],
            'camera_position': [-0.0125588,-0.00559266,0.0117482,
                                0.00438669,0.00487109,0.00135715,
                                0.334487,0.333581,0.881387],
            'color_map': 'seismic',
        }
    },
    'multi_physics/thermo_elasticity_ess.py': {
        '': {
            'fields': ['T:wu:f1e3:p0', '1:vw:p0'],
            'camera': [-45, 75, 1],
        },
    },
    'multi_physics/thermo_elasticity.py': {
        '': {'camera': [225, 75, 0.9]}
    },
    'quantum/boron.py': {
        '': {'fields': ['Psi000:p0', 'Psi001:p1', 'Psi002:p2', 'Psi003:p3'],
             'grid_vector1': [1.2, 0, 0],
             'max_plots': 4},
    },
    'quantum/hydrogen.py': {
        '': {'fields': ['Psi000:p0', 'Psi001:p1', 'Psi002:p2', 'Psi003:p3'],
             'grid_vector1': [1.2, 0, 0],
             'max_plots': 4},
    },
    'quantum/oscillator.py': {
        '': {'fields': ['Psi000:p0', 'Psi001:p1', 'Psi002:p2', 'Psi003:p3'],
             'grid_vector1': [1.2, 0, 0],
             'max_plots': 4},
    },
    'quantum/well.py': {
        '': {'fields': ['Psi000:p0', 'Psi001:p1', 'Psi002:p2', 'Psi003:p3'],
             'grid_vector1': [1.2, 0, 0],
             'max_plots': 4},
    },
}


def resview_plot(filename, filename_out, options):
    pv.set_plot_theme("document")
    plotter = pv.Plotter(off_screen=True)

    plotter = pv_plot([filename], options, plotter=plotter, use_cache=False)
    if options.axes_visibility:
        plotter.add_axes(**dict(options.axes_options))

    if options.view_2d:
        plotter.view_xy()
        plotter.show(screenshot=filename_out, window_size=(800, 600))
    else:
        if options.camera_position is not None:
            cpos = nm.array(options.camera_position)
            cpos = cpos.reshape((3, 3))
        elif options.camera:
            zoom = options.camera[2] if len(options.camera) > 2 else 1.
            cpos = get_camera_position(plotter.bounds,
                                       options.camera[0], options.camera[1],
                                       zoom=zoom)
        else:
            cpos = None

        plotter.show(cpos=cpos, screenshot=filename_out, window_size=(800, 600))


def _omit(filename, omits, omit_dirs):
    omit = False

    base = os.path.basename(filename)

    if base in omits:
        omit = True

    for omit_dir in omit_dirs:
        if omit_dir(filename) is not None:
            omit = True
            break

    return omit


def ebase2fbase(ebase):
    return os.path.splitext(ebase)[0].replace(os.path.sep, '-')


def _get_fig_filenames(ebase, images_dir):
    fig_base = ebase2fbase(ebase)

    if ebase in custom:
        suffixes = sorted(custom[ebase].keys())
        for suffix in suffixes:
            yield os.path.join(images_dir, fig_base + suffix + '.png')
    else:
        yield os.path.join(images_dir, fig_base + '.png')


def _get_fig_filename(ebase, images_dir, suffix):
    fig_base = ebase2fbase(ebase)

    return os.path.join(images_dir, fig_base + suffix + '.png')


def _make_sphinx_path(path, relative=False):
    if relative:
        aux = path.replace(sfepy.data_dir, '')
        prefix = ('..' + os.path.sep) * aux.count(os.path.sep)
        sphinx_path = prefix[:-1] + aux

    else:
        sphinx_path = path.replace(sfepy.data_dir, '/..')

    return sphinx_path


def apply_view_options(views, default):
    out = {}
    for kview, view in views.items():
        ov = default.copy()
        for k, v in view.items():
            if k == 'fields':
                fv = [k.split(':') for k in v]
                fv = [(k[0], ':'.join(k[1:])) for k in fv]
                setattr(ov, k, fv)
            else:
                setattr(ov, k, v)

        out[kview] = ov

    return out


def generate_images(images_dir, examples_dir):
    """
    Generate images from results of running examples found in
    `examples_dir` directory.

    The generated images are stored to `images_dir`,
    """
    from sfepy.applications import solve_pde
    from sfepy.solvers.ts_solvers import StationarySolver

    prefix = output.prefix

    output_dir = tempfile.mkdtemp()
    trunk = os.path.join(output_dir, 'result')
    options = Struct(output_filename_trunk=trunk,
                     output_format='vtk',
                     save_ebc=False,
                     save_ebc_nodes=False,
                     save_regions=False,
                     save_regions_as_groups=False,
                     solve_not=False)

    view_options = Struct(step=0,
                          fields=[], fields_map=[],
                          outline=False,
                          isosurfaces=0,
                          show_edges=False,
                          warp=None,
                          factor=1.,
                          opacity=1.,
                          color_map='viridis',
                          axes_options=[],
                          axes_visibility=False,
                          grid_vector1=None,
                          grid_vector2=None,
                          max_plots=3,
                          show_labels=False,
                          label_position=[-1, -1, 0, 0.2],
                          scalar_bar_size=[0.15, 0.06],
                          scalar_bar_position=[0.04, 0.92, 0, -1.5],
                          show_scalar_bars=True,
                          camera=[225, 75, 1],
                          camera_position=None,
                          view_2d=False,
                          force_view_3d=False)

    ensure_path(images_dir + os.path.sep)

    for ex_filename in locate_files('*.py', examples_dir):
        if _omit(ex_filename, omits + omit_images, omit_dirs):
            continue

        output.level = 0
        output.prefix = prefix
        ebase = ex_filename.replace(examples_dir, '')[1:]
        output('trying "%s"...' % ebase)

        try:
            problem, state = solve_pde(ex_filename, options=options)

        except KeyboardInterrupt:
            raise

        except:
            problem = None
            output('***** failed! *****')

        if problem is not None:
            if ebase in custom:
                views = apply_view_options(custom[ebase], view_options)
            else:
                views = {'': view_options.copy()}

            try:
                tsolver = problem.get_solver()

            except ValueError:
                suffix = None

            else:
                if isinstance(tsolver, StationarySolver):
                    suffix = None

                else:
                    suffix = tsolver.ts.suffix % (tsolver.ts.n_step - 1)

            filename = problem.get_output_name(suffix=suffix)
            dim = problem.get_dim()
            for suffix, kwargs in six.iteritems(views):
                if dim in (1, 2) and not kwargs.force_view_3d:
                    kwargs.view_2d = True
                    kwargs.scalar_bar_position = [0.04, 0.92, 1.7, 0]
                    if kwargs.grid_vector1 is None:
                        kwargs.grid_vector1 = [1.2, 0, 0]

                    if kwargs.grid_vector2 is None:
                        kwargs.grid_vector2 = [0, -1.2, 0]

                fig_filename = _get_fig_filename(ebase, images_dir, suffix)

                fname = edit_filename(filename, suffix=suffix)
                output('displaying results from "%s"' % fname)
                disp_name = fig_filename.replace(sfepy.data_dir, '')
                output('to "%s"...' % disp_name.lstrip(os.path.sep))

                resview_plot(fname, fig_filename, kwargs)

                output('...done')

            remove_files(output_dir)

        output('...done')


def generate_thumbnails(thumbnails_dir, images_dir, scale=0.3):
    """
    Generate thumbnails into `thumbnails_dir` corresponding to images in
    `images_dir`.
    """
    ensure_path(thumbnails_dir + os.path.sep)

    output('generating thumbnails...')
    filenames = glob.glob(os.path.join(images_dir, '*.png'))
    for fig_filename in filenames:
        ebase = fig_filename.replace(sfepy.data_dir, '').lstrip(os.path.sep)
        output('"%s"' % ebase)

        base = os.path.basename(fig_filename)
        thumb_filename = os.path.join(thumbnails_dir, base)

        image.thumbnail(fig_filename, thumb_filename, scale=scale)

    output('...done')


_index = """\
.. _%s-index:

%s
%s

.. toctree::
   :maxdepth: 2

"""

_image = '.. image:: %s'

_include = """\
.. _%s:

%s
%s

**Description**

%s

%s

:download:`source code <%s>`

.. literalinclude:: %s

"""


def generate_rst_files(rst_dir, examples_dir, images_dir):
    """
    Generate Sphinx rst files for examples in `examples_dir` with images
    in `images_dir` and put them into `rst_dir`.

    Returns
    -------
    dir_map : dict
        The directory mapping of examples and corresponding rst files.
    """
    ensure_path(rst_dir + os.path.sep)

    output('generating rst files...')

    dir_map = {}
    for ex_filename in locate_files('*.py', examples_dir):
        if _omit(ex_filename, omits, omit_dirs):
            continue

        ebase = ex_filename.replace(examples_dir, '')[1:]
        base_dir = os.path.dirname(ebase)
        rst_filename = ebase2fbase(ebase) + '.rst'
        dir_map.setdefault(base_dir, []).append((ex_filename, rst_filename))

    for dirname, filenames in six.iteritems(dir_map):
        filenames = sorted(filenames, key=lambda a: a[1])
        dir_map[dirname] = filenames

    # Main index.
    mfd = open(os.path.join(rst_dir, 'index.rst'), 'w')
    mfd.write(_index % ('examples', 'Examples', '=' * 8))

    for dirname, filenames in ordered_iteritems(dir_map):
        # Subdirectory index.
        ifd = open(os.path.join(rst_dir, '%s-index.rst' % dirname), 'w')
        ifd.write(_index % (dirname + '-examples',
                            dirname, '=' * len(dirname)))

        for ex_filename, rst_filename in filenames:
            full_rst_filename = os.path.join(rst_dir, rst_filename)
            output('"%s"' % rst_filename)
            ebase = ex_filename.replace(examples_dir, '')[1:]

            rst_base = os.path.splitext(rst_filename)[0]

            rst_ex_filename = _make_sphinx_path(ex_filename)
            docstring = get_default(import_file(ex_filename).__doc__,
                                    'missing description!')

            ifd.write('   %s <%s>\n' % (os.path.basename(ebase), rst_base))
            fig_include = ''
            for fig_filename in _get_fig_filenames(ebase, images_dir):
                rst_fig_filename = _make_sphinx_path(fig_filename)
                if os.path.exists(fig_filename):
                    fig_include += _image % rst_fig_filename + '\n'
                else:
                    output('   warning: figure "%s" not found' % fig_filename)

            # Example rst file.
            fd = open(full_rst_filename, 'w', encoding='utf-8')
            fd.write(_include % (rst_base, ebase, '=' * len(ebase),
                                 docstring,
                                 fig_include,
                                 rst_ex_filename, rst_ex_filename))
            fd.close()

        ifd.close()

        mfd.write('   %s-index\n' % dirname)

    mfd.close()

    output('...done')

    return dir_map


_rst_empty_item = """\
        - ..
"""

_rst_item = """\
      %s - .. figure:: %s
             :target: %s

             :ref:`%s <%s>`
"""

_gallery_table = """\
   .. list-table::
      :align: center
      :class: gallery
"""

_gallery_head = """\
.. only:: html

   .. _gallery-index:

   Gallery
   =======
"""


def generate_gallery(examples_dir, output_filename, doc_dir,
                     rst_dir, thumbnails_dir, dir_map, n_col=3):
    """
    Generate the gallery rst file with thumbnail images and links to
    examples.

    Parameters
    ----------
    output_filename : str
        The output rst file name.
    doc_dir : str
        The top level directory of gallery files.
    rst_dir : str
        The full path to rst files of examples within `doc_dir`.
    thumbnails_dir : str
        The full path to thumbnail images within `doc_dir`.
    dir_map : dict
        The directory mapping returned by `generate_rst_files()`
    n_col : int
        The number of columns in the gallery table.
    """
    output('generating %s...' % output_filename)

    lines = [_gallery_head]

    for dirname, filenames in ordered_iteritems(dir_map):
        title = ['   %s' % dirname.title().replace('_', ' '),
                 '   ' + len(dirname) * '^' + '',
                 _gallery_table]

        llines = []
        icol = 0
        for ex_filename, rst_filename in filenames:
            ebase = ex_filename.replace(examples_dir, '')[1:]
            link = os.path.splitext(rst_filename)[0]

            thumbnail_filename = next(_get_fig_filenames(ebase,
                                                         thumbnails_dir))
            if not os.path.isfile(thumbnail_filename):
                # Skip examples with no image (= failed examples).
                output('warning: figure "%s" not found' % thumbnail_filename)
                continue

            thumbnail_name = thumbnail_filename.replace(doc_dir, '..')
            path_to_file = os.path.join(examples_dir, ebase)
            docstring = get_default(import_file(path_to_file).__doc__,
                                    'missing description!')
            docstring = docstring.replace('e.g.', 'eg:')
            docstring = docstring.split('.')
            label = docstring[0].strip()
            label = label.replace('\n', ' ')
            label = label.replace('  ', ' ')

            llines.append(_rst_item % (' ' if icol else '*', thumbnail_name,
                                       link + '.html', label, link))
            icol = (icol + 1) % n_col

        if icol > 0:
            for j in range(icol, n_col):
                llines.append(_rst_empty_item)

        if len(llines) > 0:
            lines += title + llines
        else:
            output('warning: no figures in "%s"' % dirname)

    fd = open(output_filename, 'wt')
    fd.write('\n'.join(lines))
    fd.close()

    output('...done')


helps = {
    'doc_dir': 'top level directory of gallery files',
    'no_images': 'do not (re)generate images and thumbnails',
    'output_filename': 'output file name [default: %(default)s]',
}


def main():
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s')
    parser.add_argument('-d', '--doc-dir', metavar='doc_dir',
                        action='store', dest='doc_dir',
                        default='doc', help=helps['doc_dir'])
    parser.add_argument('-n', '--no-images',
                        action='store_true', dest='no_images',
                        default=False, help=helps['no_images'])
    parser.add_argument('-o', '--output', metavar='output_filename',
                        action='store', dest='output_filename',
                        default='gallery.rst',
                        help=helps['output_filename'])
    options = parser.parse_args()

    rst_dir = 'examples'
    doc_dir = os.path.realpath(options.doc_dir)
    examples_dir = os.path.realpath('sfepy/examples')
    full_rst_dir = os.path.join(doc_dir, rst_dir)
    images_dir = os.path.join(doc_dir, 'images/gallery')
    thumbnails_dir = os.path.join(images_dir, 'thumbnails')

    output_filename = os.path.join(full_rst_dir, options.output_filename)

    if not options.no_images:
        generate_images(images_dir, examples_dir)
        generate_thumbnails(thumbnails_dir, images_dir)

    dir_map = generate_rst_files(full_rst_dir, examples_dir, images_dir)

    generate_gallery(examples_dir, output_filename, doc_dir,
                     rst_dir, thumbnails_dir, dir_map)


if __name__ == '__main__':
    main()

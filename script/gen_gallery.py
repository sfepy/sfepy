#!/usr/bin/env python
"""
Generate the images and rst files for gallery of SfePy examples.

The following steps need to be made to regenerate the documentation with the
updated example files:

1. remove doc/examples/*::

   $ rm -rf doc/examples/*

2. generate the files:

   $ ./script/gen_gallery.py

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

import matplotlib.image as image

import sfepy
from sfepy.base.base import (get_default, ordered_iteritems,
                             import_file, output, Struct)
from sfepy.base.ioutils import (ensure_path, locate_files, remove_files,
                                edit_filename)
from sfepy.postprocess.domain_specific import DomainSpecificPlot

omits = [
    'vibro_acoustic3d_mid.py',
    'its2D_5.py',
    'linear_elastic_probes.py',
    '__init__.py',
]

omit_dirs = [
    re.compile('.*output.*/').match,
]

custom = {
    'acoustics/acoustics3d.py' : {
        '_p_1' : {
            'is_scalar_bar' : True,
            'view' : (44, 57, 0.24, [-0.004, -0.007, 0.09]),
            'roll' : 0,
        },
        '_p_2' : {
            'is_scalar_bar' : True,
            'view' : (-99, 120, 0.4, [0.0, 0.0, 0.07]),
            'roll' : 141,
        },
    },
    'acoustics/vibro_acoustic3d.py' : {
        '_p1' : {
            'is_scalar_bar' : True,
            'view' : (45.0, 54.7, 1.47, [0.0, 0.0, 0.05]),
            'roll' : -120,
        },
        '_p2' : {
            'is_scalar_bar' : True,
            'view' : (45.0, 54.7, 1.47, [0.0, 0.0, 0.15]),
            'roll' : -120,
        },
        '_w' : {
            'is_scalar_bar' : True,
            'view' : (0.0, 0.0, 0.86, [0.0, 0.0, 0.1]),
            'roll' : 0,
        },
        '_g0' : {
            'is_scalar_bar' : True,
            'view' : (0.0, 0.0, 0.86, [0.0, 0.0, 0.1]),
            'roll' : 0,
        },
    },
    'diffusion/laplace_1d.py' : {
        '' : {
            'is_scalar_bar' : True,
            'is_wireframe' : True,
            'domain_specific' : {
                't' : DomainSpecificPlot('plot_warp_scalar',
                                          ['rel_scaling=1']),
            },
            'view' : (-90, 90, 1.5, [0,  0, 0]),
            'roll' : 0,
            'opacity' : {'wireframe' : 0.3},
        },
    },
    'diffusion/laplace_coupling_lcbcs.py' : {
        '' : {
            'is_scalar_bar' : True,
            'is_wireframe' : True,
            'domain_specific' : {
                'u1' : DomainSpecificPlot('plot_warp_scalar',
                                          ['rel_scaling=1']),
                'u2' : DomainSpecificPlot('plot_warp_scalar',
                                          ['rel_scaling=1']),
            },
            'view' : (-82, 50, 3.6, [-0.43, -0.55, 0.4]),
            'roll' : -23,
            'opacity' : {'wireframe' : 0.3},
        },
    },
    'diffusion/poisson_iga.py' : {
        '' : {
            'is_scalar_bar' : True,
            'is_wireframe' : True,
            'domain_specific' : {
                't' : DomainSpecificPlot('plot_warp_scalar',
                                         ['rel_scaling=1']),
            },
            'view' : (55, 39, 6.6, [-0.35, -0.29, 0.35]),
            'roll' : 15,
            'opacity' : {'wireframe' : 0.3},
        },
    },
    'diffusion/sinbc.py' : {
        '_t' : {
            'is_scalar_bar' : True,
            'is_wireframe' : True,
            'domain_specific' : {
                't' : DomainSpecificPlot('plot_warp_scalar',
                                         ['rel_scaling=1']),
            },
            'view' : (-170, 30, 4.7, [0.34, 0.23, -0.26]),
            'roll' : 71,
            'opacity' : {'wireframe' : 0.3},
        },
        '_grad' : {
            'is_scalar_bar' : True,
            'opacity' : {'surface' : 0.3},
            'view' : (-170, 30, 4.7, [0.34, 0.23, -0.26]),
            'roll' : 71,
        },
    },
    'linear_elasticity/elastic_contact_planes.py' : {
        '' : {
            'is_scalar_bar' : True,
            'is_wireframe' : True,
            'domain_specific' : {
                'u' : DomainSpecificPlot('plot_displacements',
                                         ['rel_scaling=1']),
            },
            'view' : (-82, 47, 3.4, [-0.5, -0.24, -0.2]),
            'roll' : -8.4,
            'opacity' : {'wireframe' : 0.3},
        },
    },
    'linear_elasticity/elastic_contact_sphere.py' : {
        '' : {
            'is_scalar_bar' : True,
            'is_wireframe' : True,
            'domain_specific' : {
                'u' : DomainSpecificPlot('plot_displacements',
                                         ['rel_scaling=1']),
            },
            'view' : (-82, 47, 3.4, [-0.5, -0.24, -0.2]),
            'roll' : -8.4,
            'opacity' : {'wireframe' : 0.3},
        },
    },
    'linear_elasticity/elastic_shifted_periodic.py' : {
        '' : {
            'is_scalar_bar' : True,
            'is_wireframe' : True,
            'only_names' : ['u'],
            'domain_specific' : {
                'u' : DomainSpecificPlot('plot_displacements',
                                         ['rel_scaling=1',
                                          'color_kind="scalars"',
                                          'color_name="von_mises_stress"']),
            },
            'view' : (142, 39, 16, [-4.7, -2.1, -1.9]),
            'roll' : 8.4,
            'opacity' : {'wireframe' : 0.3},
        },
    },
    'linear_elasticity/linear_elastic_iga.py' : {
        '' : {
            'is_scalar_bar' : True,
            'is_wireframe' : True,
            'domain_specific' : {
                'u' : DomainSpecificPlot('plot_displacements',
                                         ['rel_scaling=1']),
            },
            'view' : (-37, 51, 1.5, [-0.28, -0.29, 0.0]),
            'roll' : -51.5,
            'opacity' : {'wireframe' : 0.2},
        },
    },
    'linear_elasticity/shell10x_cantilever.py' : {
        '' : {
            'is_scalar_bar' : True,
            'is_wireframe' : True,
            'domain_specific' : {
                'u_disp' : DomainSpecificPlot('plot_displacements',
                                              ['rel_scaling=1']),
            },
            'view' : (-45, 81, 0.59, [-0.075,  0.023,  0.093]),
            'roll' : -75.0,
            'opacity' : {'wireframe' : 0.5},
        },
    },
    'navier_stokes/stokes_slip_bc.py' : {
        '' : {
            'is_scalar_bar' : True,
            'view' : (-63, 52, 5.2, [-1.5, -0.65, 0.12]),
            'roll' : -32,
            'resolution' : (800, 600),
            'layout' : 'col',
            'rel_scaling' : 0.1,
        },
    },
    'navier_stokes/stokes_slip_bc_penalty.py' : {
        '' : {
            'is_scalar_bar' : True,
            'view' : (-63, 52, 5.2, [-1.5, -0.65, 0.12]),
            'roll' : -32,
            'resolution' : (800, 600),
            'layout' : 'col',
            'rel_scaling' : 0.1,
        },
    },
    'multi_physics/thermo_elasticity_ess.py' : {
        '' : {
            'is_scalar_bar' : True,
            'is_wireframe' : True,
            'only_names' : ['u'],
            'domain_specific' : {
                'u' : DomainSpecificPlot('plot_displacements',
                                         ['rel_scaling=1000',
                                          'color_kind="scalars"',
                                          'color_name="T"']),
            },
            'view' : (-51, 71, 12.9, [-2.3, -2.4, -0.2]),
            'roll' : -65,
            'opacity' : {'wireframe' : 0.3},
        },
    },
    'quantum/boron.py' : {
        '' : {
            'is_scalar_bar' : False,
            'only_names' : ['Psi000', 'Psi001', 'Psi002'],
            'view' : (0.0, 0.0, 200.0, [-25., -25.,   0.]),
            'roll' : 0.0,
        },
    },
    'quantum/hydrogen.py' : {
        '' : {
            'is_scalar_bar' : False,
            'only_names' : ['Psi000', 'Psi001', 'Psi002'],
            'view' : (0.0, 0.0, 200.0, [-25., -25.,   0.]),
            'roll' : 0.0,
        },
    },
    'quantum/oscillator.py' : {
        '' : {
            'is_scalar_bar' : False,
            'only_names' : ['Psi000', 'Psi001', 'Psi002'],
            'view' : (0.0, 0.0, 200.0, [-25., -25.,   0.]),
            'roll' : 0.0,
        },
    },
    'quantum/well.py' : {
        '' : {
            'is_scalar_bar' : False,
            'only_names' : ['Psi000', 'Psi001', 'Psi002'],
            'view' : (0.0, 0.0, 200.0, [-25., -25.,   0.]),
            'roll' : 0.0,
        },
    },
}

def _omit(filename):
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

def generate_images(images_dir, examples_dir):
    """
    Generate images from results of running examples found in
    `examples_dir` directory.

    The generated images are stored to `images_dir`,
    """
    from sfepy.applications import solve_pde
    from sfepy.postprocess.viewer import Viewer
    from sfepy.postprocess.utils import mlab
    from sfepy.solvers.ts_solvers import StationarySolver

    prefix = output.prefix

    output_dir = tempfile.mkdtemp()
    trunk = os.path.join(output_dir, 'result')
    options = Struct(output_filename_trunk=trunk,
                     output_format='vtk',
                     save_ebc=False,
                     save_ebc_nodes=False,
                     save_regions=False,
                     save_field_meshes=False,
                     save_regions_as_groups=False,
                     solve_not=False)
    default_views = {'' : {'is_scalar_bar' : True}}

    ensure_path(images_dir + os.path.sep)

    view = Viewer('', offscreen=False)

    for ex_filename in locate_files('*.py', examples_dir):
        if _omit(ex_filename): continue

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
                views = custom[ebase]

            else:
                views = default_views

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
            for suffix, kwargs in six.iteritems(views):
                fig_filename = _get_fig_filename(ebase, images_dir, suffix)

                fname = edit_filename(filename, suffix=suffix)
                output('displaying results from "%s"' % fname)
                disp_name = fig_filename.replace(sfepy.data_dir, '')
                output('to "%s"...' % disp_name.lstrip(os.path.sep))

                view.filename = fname
                view(scene=view.scene, show=False, **kwargs)
                view.save_image(fig_filename)
                mlab.clf()

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
        if _omit(ex_filename): continue

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

            rst_base = rst_filename.replace('.rst', '')

            rst_ex_filename = _make_sphinx_path(ex_filename)
            docstring = get_default(import_file(ex_filename).__doc__,
                                    'missing description!')

            ifd.write('    %s <%s>\n' % (os.path.basename(ebase),
                                         rst_filename.replace('.rst', '')))
            fig_include = ''
            for fig_filename in _get_fig_filenames(ebase, images_dir):
                rst_fig_filename = _make_sphinx_path(fig_filename)
                if os.path.exists(fig_filename):
                    fig_include += _image % rst_fig_filename + '\n'
                else:
                    output('   warning: figure "%s" not found' % fig_filename)
 

            # Example rst file.
            fd = open(full_rst_filename, 'w')
            fd.write(_include % (rst_base, ebase, '=' * len(ebase),
                                 docstring,
                                 fig_include,
                                 rst_ex_filename, rst_ex_filename))
            fd.close()

        ifd.close()

        mfd.write('    %s-index\n' % dirname)

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
        title= ['%s' % dirname.title().replace('_', ' '),
                len(dirname) * '^' + '',
                _gallery_table]

        llines = []
        icol = 0
        for ex_filename, rst_filename in filenames:
            ebase = ex_filename.replace(examples_dir, '')[1:]
            link = rst_filename.replace('.rst', '')

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
    examples_dir = os.path.realpath('examples')
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

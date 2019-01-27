#!/usr/bin/env python
"""
Generate the images and rst files for gallery of SfePy examples.

The following steps need to be made to regenerate the documentation with the
updated example files:

1. Generate the files:

   - for sfepy.org deployment::

     $ ./script/gen_gallery.py -l ../doc-devel

   - for local test build run from ./::

     $ ./script/gen_gallery.py -l doc/_build/html/

2. remove doc/examples/::

   $ rm -rf doc/examples/

3. copy gallery/examples/ to doc/::

   $ cp -a gallery/examples/ doc/

4. regenerate the documentation::

   $ python setup.py htmldocs

Additional steps for sfepy.org deployment:

- copy doc/_build/html/ to <sfepy.org>/doc-devel/
- copy gallery/index.html and gallery/images/ to <sfepy.org>/
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

def _get_fig_filenames(ebase, images_dir):
    fig_base = os.path.splitext(ebase)[0].replace(os.path.sep, '-')

    yield fig_base

    if ebase in custom:
        suffixes = sorted(custom[ebase].keys())
        for suffix in suffixes:
            fig_filename = os.path.join(images_dir, fig_base + suffix + '.png')
            yield fig_filename

    else:
        fig_filename = os.path.join(images_dir, fig_base + '.png')
        yield fig_filename

def _get_fig_filename(ebase, images_dir, suffix):
    fig_base = os.path.splitext(ebase)[0].replace(os.path.sep, '-')
    fig_filename = os.path.join(images_dir, fig_base + suffix + '.png')

    return fig_filename

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
.. _%s-gallery-examples-index:

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

        rst_filename = os.path.basename(ex_filename).replace('.py', '.rst')

        dir_map.setdefault(base_dir, []).append((ex_filename, rst_filename))

    for dirname, filenames in six.iteritems(dir_map):
        filenames = sorted(filenames, key=lambda a: a[1])
        dir_map[dirname ] = filenames

    # Main index.
    mfd = open(os.path.join(rst_dir, 'index.rst'), 'w')
    mfd.write(_index % ('sfepy', 'Examples', '=' * 8))

    for dirname, filenames in ordered_iteritems(dir_map):
        full_dirname = os.path.join(rst_dir, dirname)
        ensure_path(full_dirname + os.path.sep)

        # Subdirectory index.
        ifd = open(os.path.join(full_dirname, 'index.rst'), 'w')
        ifd.write(_index % (dirname, dirname, '=' * len(dirname)))

        for ex_filename, rst_filename in filenames:
            full_rst_filename = os.path.join(full_dirname, rst_filename)
            output('"%s"' % full_rst_filename.replace(rst_dir, '')[1:])
            rst_filename_ns = rst_filename.replace('.rst', '')
            ebase = ex_filename.replace(examples_dir, '')[1:]

            rst_ex_filename = _make_sphinx_path(ex_filename)
            docstring = get_default(import_file(ex_filename).__doc__,
                                    'missing description!')

            ifd.write('    %s\n' % rst_filename_ns)
            fig_include = ''
            fig_base = next(_get_fig_filenames(ebase, images_dir))
            for fig_filename in _get_fig_filenames(ebase, images_dir):
                rst_fig_filename = _make_sphinx_path(fig_filename)

                if os.path.exists(fig_filename):
                    fig_include += _image % rst_fig_filename + '\n'

            # Example rst file.
            fd = open(full_rst_filename, 'w')
            fd.write(_include % (fig_base, ebase, '=' * len(ebase),
                                 docstring,
                                 fig_include,
                                 rst_ex_filename, rst_ex_filename))
            fd.close()

        ifd.close()

        mfd.write('    %s/index\n' % dirname)

    mfd.close()

    output('...done')

    return dir_map

_gallery_template_file = os.path.join(sfepy.top_dir,
                                      'doc/gallery_template.html')

_link_template = """\
<div class="figure">
<a class="reference external image-reference" href="../%s">
<img alt="%s" src="%s" />
</a>
<p class="caption">
<a class="reference internal" href="../%s"><em>%s</em></a>
</p>
</div>
<div class="toctree-wrapper compound">
</div>
"""
_side_links="<li><a class='reference internal' href='#%s'>%s</a></li>"
_div_line ="""\
<div class="section" id="%s">
<h2>%s<a class="headerlink" href="\#%s" title="Permalink to this headline">
</a></h2>
%s
<div style="clear: both"></div></div>
"""
def generate_gallery_html(examples_dir, output_filename, gallery_dir,
                          rst_dir, thumbnails_dir, dir_map, link_prefix):
    """
    Generate the gallery html file with thumbnail images and links to
    examples.

    Parameters
    ----------
    output_filename : str
        The output html file name.
    gallery_dir : str
        The top level directory of gallery files.
    rst_dir : str
        The full path to rst files of examples within `gallery_dir`.
    thumbnails_dir : str
        The full path to thumbnail images within `gallery_dir`.
    dir_map : dict
        The directory mapping returned by `generate_rst_files()`
    link_prefix : str, optional
        The prefix to prepend to links to individual pages of examples.
    """
    output('generating %s...' % output_filename)

    with open(_gallery_template_file, 'r') as fd:
        gallery_template = fd.read()

    div_lines=[]
    sidebar = []
    for dirname, filenames in ordered_iteritems(dir_map):
        full_dirname = os.path.join(rst_dir, dirname)
        dirnamenew = dirname.replace("_"," ")
        sidebarline = _side_links % (dirname, dirnamenew.title())
        lines = []
        for ex_filename, rst_filename in filenames:
            full_rst_filename = os.path.join(full_dirname, rst_filename)

            ebase = full_rst_filename.replace(rst_dir, '')[1:]
            ebase = edit_filename(ebase, new_ext='.py')

            link_base = full_rst_filename.replace(gallery_dir, '')[1:]
            link = os.path.join(link_prefix,
                                os.path.splitext(link_base)[0] + '.html')

            next(_get_fig_filenames(ebase, thumbnails_dir))
            for thumbnail_filename in _get_fig_filenames(ebase,
                                                         thumbnails_dir):
                if not os.path.isfile(thumbnail_filename):
                    # Skip examples with no image (= failed examples).
                    continue

                thumbnail_name = thumbnail_filename.replace(gallery_dir,
                                                            '')[1:]
                path_to_file = os.path.join(examples_dir,ebase)
                docstring = get_default(import_file(path_to_file).__doc__,
                                        'missing description!')
                docstring = docstring.replace('e.g.', 'eg:')
                docstring = docstring.split('.')
                line = _link_template % (link,os.path.splitext(ebase)[0],
                                         thumbnail_name,link,docstring[0]+'.')
                lines.append(line)

        if(len(lines)!=0):
            div_lines.append(_div_line % (dirname, dirnamenew.title(),
                                          dirname, '\n'.join(lines)))
            sidebar.append(sidebarline)

    fd = open(output_filename, 'w')
    fd.write(gallery_template % ((link_prefix,) * 7
                                 + ('\n'.join(sidebar), '\n'.join(div_lines))))
    fd.close()

    output('...done')

helps = {
    'examples_dir' :
    'directory containing examples [default: %(default)s]',
    'images_dir' :
    'directory where to store gallery images [default: gallery/images]',
    'no_images' :
    'do not (re)generate images and thumbnails',
    'output_filename' :
    'output file name [default: %(default)s]',
    'link_prefix' :
    'prefix to be prepended to links to examples pages in gallery '
    '[default: %(default)s]',
}

def main():
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s')
    parser.add_argument('-e', '--examples-dir', metavar='directory',
                        action='store', dest='examples_dir',
                        default='examples', help=helps['examples_dir'])
    parser.add_argument('-i', '--images-dir', metavar='directory',
                        action='store', dest='images_dir',
                        default=None, help=helps['images_dir'])
    parser.add_argument('-n', '--no-images',
                        action='store_true', dest='no_images',
                        default=False, help=helps['no_images'])
    parser.add_argument('-o', '--output', metavar='output_filename',
                        action='store', dest='output_filename',
                        default='gallery/index.html',
                        help=helps['output_filename'])
    parser.add_argument('-l', '--link-prefix', metavar='prefix',
                        action='store', dest='link_prefix',
                        default='http://sfepy.org/doc-devel',
                        help=helps['link_prefix'])
    options = parser.parse_args()

    examples_dir = os.path.realpath(options.examples_dir)

    output_filename = os.path.realpath(options.output_filename)
    gallery_dir = os.path.dirname(output_filename)

    images_dir = get_default(options.images_dir,
                             os.path.join(gallery_dir, 'images'))

    thumbnails_dir = os.path.join(images_dir, 'thumbnails')
    rst_dir = os.path.join(gallery_dir, 'examples')
    if not options.no_images:
        generate_images(images_dir, examples_dir)
        generate_thumbnails(thumbnails_dir, images_dir)

    dir_map = generate_rst_files(rst_dir, examples_dir, images_dir)

    generate_gallery_html(examples_dir,output_filename, gallery_dir,
                          rst_dir, thumbnails_dir, dir_map,
                          link_prefix=options.link_prefix)

if __name__ == '__main__':
    main()

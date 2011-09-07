#!/usr/bin/env python
import sys
sys.path.append( '.' )
import os
import tempfile
from optparse import OptionParser

import sfepy
from sfepy.base.base import output, Struct
from sfepy.base.ioutils import ensure_path, locate_files, remove_files
from sfepy.applications import pde_solve
from sfepy.postprocess import Viewer
from sfepy.postprocess.utils import mlab

omits = [
    'linear_elastic_mM.py',
]

def generate_images(examples_dir, image_dir='gallery/images'):
    """
    Generate images from results of running examples found in
    `examples_dir` directory.

    The generated images are stored to `image_dir`,
    """
    prefix = output.prefix

    output_dir = tempfile.mkdtemp()
    trunk = os.path.join(output_dir, 'result')
    options = Struct(output_filename_trunk=trunk,
                     output_format='vtk',
                     save_ebc=False,
                     save_regions=False,
                     save_field_meshes=False,
                     save_regions_as_groups=False,
                     solve_not=False)

    cwd = os.path.dirname(os.path.realpath(examples_dir))
    gallery_dir = os.path.join(cwd, image_dir + os.path.sep)
    ensure_path(gallery_dir)

    view = Viewer('',
                  output_dir=output_dir,
                  offscreen=False)

    for ex_filename in locate_files('*.py', examples_dir):
        base = os.path.basename(ex_filename)
        if base in omits: continue

        output.level = 0
        output.prefix = prefix
        ebase = ex_filename.replace(sfepy.data_dir, '')[1:]
        output('trying "%s"...' % ebase)

        try:
            problem, state = pde_solve(ex_filename, options=options)

        except KeyboardInterrupt:
            raise

        except:
            problem = None
            output('***** failed! *****')

        if problem is not None:
            fig_base = ebase[:-3].replace(os.path.sep, '#') + '.png'
            fig_filename = os.path.join(gallery_dir, fig_base)

            if problem.ts_conf is None:
                filename = trunk + '.vtk'

            else:
                suffix = problem.ts.suffix % problem.ts.step
                filename = problem.get_output_name(suffix=suffix)

            output('displaying results from "%s"' % filename)
            output('to "%s"...'
                   % fig_filename.replace(sfepy.data_dir, '')[1:])

            view.filename = filename
            view(scene=view.scene, show=False, is_scalar_bar=True,
                 fig_filename=fig_filename)
            mlab.clf()

            output('...done')

            remove_files(output_dir)

        output('...done')

usage = """%prog [options]

Generate the images for gallery of SfePy examples.
"""
help = {
    'output_filename' :
    'output file name [default: %default]',
    'examples_dir' :
    'directory containing examples [default: %default]',
}

def main():
    parser = OptionParser(usage=usage, version='%prog')
    parser.add_option('-e', '--examples-dir', metavar='directory',
                      action='store', dest='examples_dir',
                      default='examples', help=help['examples_dir'])
    (options, args) = parser.parse_args()

    generate_images(options.examples_dir)

if __name__ == '__main__':
    main()

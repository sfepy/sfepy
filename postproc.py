#!/usr/bin/env python
from optparse import OptionParser
import os
import glob

import sfepy
from sfepy.base.base import pause, output
from sfepy.postprocess import Viewer, get_data_ranges, create_file_source
from sfepy.solvers.ts import get_print_info

usage = """%prog [options] filename

This is a script for quick Mayavi-based visualizations of finite element
computations results.

Examples
--------
  The examples assume that runTests.py has been run successfully and the
  resulting data files are present.

  - view data in output-tests/test_navier_stokes.vtk

    $ python postproc.py output-tests/test_navier_stokes.vtk
    $ python postproc.py output-tests/test_navier_stokes.vtk --3d

  - create animation (forces offscreen rendering) from
    output-tests/test_time_poisson.*.vtk
    
    $ python postproc.py output-tests/test_time_poisson.*.vtk -a mov

  - create animation (forces offscreen rendering) from
    output-tests/test_hyperelastic.*.vtk

    The range specification for the displacements 'u' is required, as
    output-tests/test_hyperelastic.00.vtk contains only zero
    displacements which leads to invisible glyph size.
    
    $ python postproc.py output-tests/test_hyperelastic.*.vtk \
                         --ranges=u,0,0.02 -a mov 

  - same as above, but slower frame rate

    $ python postproc.py output-tests/test_hyperelastic.*.vtk \
                         --ranges=u,0,0.02 -a mov --ffmpeg-options="-r 2 -sameq"

"""

help = {
    'list_ranges' :
    'do not plot, only list names and ranges of all data',
    'no_show' :
    'do not call mlab.show()',
    'is_3d' :
    '3d plot mode',
    'view' :
    'camera view angles [default: if --3d is True: "45,45", else: "0,0"]',
    'roll' :
    'camera roll angle [default: %default]',
    'layout' :
    'layout for multi-field plots, one of: rowcol, colrow, row, col' \
    ' [default: %default]',
    'scalar_mode' :
    'mode for plotting scalars with --3d, one of: cut_plane, iso_surface,'\
    ' both [default: %default]',
    'rel_scaling' :
    'relative scaling of glyphs (vector field visualization)' \
    ' [default: %default]',
    'clamping' :
    'glyph clamping mode',
    'ranges' :
    'force data ranges [default: automatic from data]',
    'is_scalar_bar' :
    'show scalar bar for each data',
    'rel_text_width' :
    'relative text annotation width [default: %default]',
    'filename' :
    'view image file name [default: %default]',
    'anim_file_type' :
    'if set to a ffmpeg-supported format (e.g. mov, avi, mpg), ffmpeg is' \
    ' installed and results of multiple time steps are given, an animation is' \
    ' created in the same directory as the view images',
    'ffmpeg_options' :
    'ffmpeg animation encoding options (enclose in "") [default: %default]',
    'resolution' :
    'image resolution in NxN format [default: shorter axis: 600;'\
    ' depends on layout: for rowcol it is 800x600]',
    'all' :
    'draw all data (normally, node_groups and mat_id are omitted)',
    'only_names' :
    'draw only named data',
    'anti_aliasing' :
    'value of anti-aliasing [default: mayavi2 default]',
}

def parse_view(option, opt, value, parser):
    val = tuple(float(ii) for ii in value.split(','))
    setattr(parser.values, option.dest, val)

def parse_resolution(option, opt, value, parser):
    if value is not None:
        print value
        setattr(parser.values, option.dest,
                tuple([int(r) for r in value.split('x')]))

def parse_ranges(option, opt, value, parser):
    if value is not None:
        print value
        ranges = {}
        for rng in value.split(':'):
            aux = rng.split(',')
            ranges[aux[0]] = (float(aux[1]), float(aux[2]))
        setattr(parser.values, option.dest, ranges)

def view_single_file(filename, filter_names, options, view=None):
    if view is None:
        view = Viewer(filename, offscreen=not options.show)

        if options.list_ranges:
            file_source = create_file_source(filename)
            file_source.set_step(100)
            get_data_ranges(file_source())

        else:
            if options.only_names is not None:
                options.only_names = options.only_names.split(',')

            view(show=options.show, is_3d=options.is_3d, view=options.view,
                 roll=options.roll, layout=options.layout,
                 scalar_mode=options.scalar_mode,
                 rel_scaling=options.rel_scaling,
                 clamping=options.clamping, ranges=options.ranges,
                 is_scalar_bar=options.is_scalar_bar,
                 rel_text_width=options.rel_text_width,
                 fig_filename=options.filename, resolution=options.resolution,
                 filter_names=filter_names, only_names=options.only_names,
                 anti_aliasing=options.anti_aliasing)
    else:
        view.set_source_filename(filename)
        view.save_image(options.filename)

    return view

                
def main():
    parser = OptionParser(usage=usage, version="%prog " + sfepy.__version__)
    parser.add_option("-l", "--list-ranges",
                      action="store_true", dest="list_ranges",
                      default=False, help=help['list_ranges'])
    parser.add_option("-n", "--no-show",
                      action="store_false", dest="show",
                      default=True, help=help['no_show'])
    parser.add_option("--3d",
                      action="store_true", dest="is_3d",
                      default=False, help=help['is_3d'])
    parser.add_option("--view", type='str', metavar='angle,angle',
                      action="callback", dest="view",
                      callback=parse_view, help=help['view'])
    parser.add_option("--roll", type='float', metavar='angle',
                      action="store", dest="roll",
                      default=0.0, help=help['roll'])
    parser.add_option("--layout", metavar='layout',
                      action="store", dest="layout",
                      default='rowcol', help=help['layout'])
    parser.add_option("--scalar-mode", metavar='mode',
                      action="store", dest="scalar_mode",
                      default='iso_surface', help=help['scalar_mode'])
    parser.add_option("-s", "--scale-glyphs", type='float', metavar='scale',
                      action="store", dest="rel_scaling",
                      default=0.05, help=help['rel_scaling'])
    parser.add_option("--clamping",
                      action="store_true", dest="clamping",
                      default=False, help=help['clamping'])
    parser.add_option("--ranges", type='str',
                      metavar='name1,min1,max1:name2,min2,max2:...',
                      action="callback", dest="ranges",
                      callback=parse_ranges, help=help['ranges'])
    parser.add_option("-b", "--scalar-bar",
                      action="store_true", dest="is_scalar_bar",
                      default=False, help=help['is_scalar_bar'])
    parser.add_option("--rel-text-width", type='float', metavar='width',
                      action="store", dest="rel_text_width",
                      default=0.02, help=help['rel_text_width'])
    parser.add_option("-o", "--output", metavar='filename',
                      action="store", dest="filename",
                      default='view.png', help=help['filename'])
    parser.add_option("-a", "--animation", metavar='<ffmpeg-supported format>',
                      action="store", dest="anim_file_type",
                      default=None, help=help['anim_file_type'])
    parser.add_option("", "--ffmpeg-options", metavar='"<ffmpeg options>"',
                      action="store", dest="ffmpeg_options",
                      default='-r 10 -sameq',
                      help=help['ffmpeg_options'])
    parser.add_option("-r", "--resolution", type='str', metavar='resolution',
                      action="callback", dest="resolution",
                      callback=parse_resolution, help=help['resolution'])
    parser.add_option("--all",
                      action="store_true", dest="all",
                      default=False, help=help['all'])
    parser.add_option("--only-names", metavar='list of names',
                      action="store", dest="only_names",
                      default=None, help=help['only_names'])
    parser.add_option("--anti-aliasing", type='int', metavar='value',
                      action="store", dest="anti_aliasing",
                      default=None, help=help['anti_aliasing'])
    options, args = parser.parse_args()

    if len(args) >= 1:
        if len(args) == 1:
            filenames = glob.glob(args[0])
            filenames.sort()
        else:
            filenames = args
    else:
        parser.print_help(),
        return

    if not options.all:
        filter_names = ['node_groups', 'mat_id']
    else:
        filter_names = []

    if len(filenames) == 1:
        view_single_file(filenames[0], filter_names, options)        

    else:
        if options.anim_file_type is not None:
            # Force the offscreen rendering when saving an animation.
            options.show = False

        fig_filename=options.filename
        base, ext = os.path.splitext(fig_filename)

        n_digit, fmt, suffix = get_print_info(len(filenames))

        view = None
        for ii, filename in enumerate(filenames):
            output('%d: %s' % (ii, filename))
            options.filename = '.'.join((base, suffix % ii, ext[1:]))
            view = view_single_file(filename, filter_names, options, view)

        if options.anim_file_type is not None:
            anim_name = '.'.join((base, options.anim_file_type))
            cmd = 'ffmpeg %s -i %s %s' % (options.ffmpeg_options,
                                          '.'.join((base, suffix, ext[1:])),
                                          anim_name)
            output('creating animation "%s"...' % anim_name)
            try:
                os.system(cmd) 
            except:
                output('...warning: animation not created, is ffmpeg installed?')
            else:
                output('...done')

if __name__ == '__main__':
    main()

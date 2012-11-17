#!/usr/bin/env python
from optparse import OptionParser
import os
import glob

import sfepy
from sfepy.base.base import assert_, get_default, output, nm
from sfepy.postprocess import Viewer, get_data_ranges, create_file_source
from sfepy.postprocess.domain_specific import DomainSpecificPlot

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
    'step' :
    'set the time step, negative indices are allowed, -1 means the last step'
    ' [default: %default]',
    'no_show' :
    'do not call mlab.show()',
    'no_offscreen' :
    'force no offscreen rendering for --no-show',
    'is_3d' :
    '3d plot mode',
    'view' :
    'camera azimuth, elevation angles, and optionally also '
    'distance and focal point coordinates (without []) as in `mlab.view()` '
    '[default: if --3d is True: "45,45", else: "0,0"]',
    'roll' :
    'camera roll angle [default: %default]',
    'fgcolor' :
    'foreground color, that is the color of all text annotation labels'
    ' (axes, orientation axes, scalar bar labels) [default: %default]',
    'bgcolor' :
    'background color [default: %default]',
    'layout' :
    'layout for multi-field plots, one of: rowcol, colrow, row, col' \
    ' [default: %default]',
    'scalar_mode' :
    'mode for plotting scalars with --3d, one of: cut_plane, iso_surface,'\
    ' both [default: %default]',
    'vector_mode' :
    'mode for plotting vectors, one of: arrows, norm, arrows_norm, warp_norm'\
    ' [default: %default]',
    'rel_scaling' :
    'relative scaling of glyphs (vector field visualization)' \
    ' [default: %default]',
    'clamping' :
    'glyph clamping mode',
    'ranges' :
    'force data ranges [default: automatic from data]',
    'is_scalar_bar' :
    'show scalar bar for each data',
    'is_wireframe' :
    'show wireframe of mesh surface for each data',
    'opacity' :
    'opacity in [0.0, 1.0]. Can be given either globally'
    ' as a single float, or per module, e.g.'
    ' "wireframe=0.1,scalar_cut_plane=0.5". Possible keywords are: wireframe,'
    ' scalar_cut_plane, vector_cut_plane, surface, iso_surface,'
    ' arrows_surface, glyphs. [default: 1.0]',
    'rel_text_width' :
    'relative text annotation width [default: %default]',
    'watch' :
    'watch the results file for changes (single file mode only)',
    'filename' :
    "view image file name [default: 'view.png']",
    'output_dir' :
    "output directory for saving view images; ignored when -o option is" \
    " given, as the directory part of the filename is taken instead" \
    " [default: '.']",
    'anim_format' :
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
    'group_names' :
    'superimpose plots of data in each group',
    'subdomains' :
    'superimpose surfaces of subdomains over each data;' \
    ' example value: mat_id,0,None,True',
    'anti_aliasing' :
    'value of anti-aliasing [default: mayavi2 default]',
    'domain_specific' :
    'domain specific drawing functions and configurations',
}

def parse_view(option, opt, value, parser):
    vals = value.split(',')
    assert_(len(vals) in [2, 3, 6])
    val = tuple(float(ii) for ii in vals)
    if len(vals) == 6:
        val = val[:3] + (list(val[3:]),)
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

def parse_opacity(option, opt, value, parser):
    try:
        opacity = float(value)
        assert_(0.0 <= opacity <= 1.0)

    except:
        opacity = {}

        for vals in value.split(','):
            key, val = vals.split('=')
            val = float(val)
            assert_(0.0 <= val <= 1.0)

            opacity[key] = val

    setattr(parser.values, option.dest, opacity)

def parse_group_names(option, opt, value, parser):
    if value is not None:
        print value
        group_names = [tuple(group.split(',')) for group in value.split(':')]
        setattr(parser.values, option.dest, group_names)

def parse_subdomains(option, opt, value, parser):
    if value is not None:
        print value
        aux = value.split(',')

        try:
            tmin = int(aux[1])

        except ValueError:
            tmin = None

        try:
            tmax = int(aux[2])

        except ValueError:
            tmax = None

        subdomains_args = {'mat_id_name' : aux[0],
                           'threshold_limits' : (tmin, tmax),
                           'single_color' : aux[3] == 'True'}
        setattr(parser.values, option.dest, subdomains_args)

def parse_domain_specific(option, opt, value, parser):
    if value is not None:
        print value
        out = {}
        confs = value.split(':')
        for conf in confs:
            aux = conf.split(',')
            var_name, fun_name = aux[:2]
            args = aux[2:]

            out[var_name] = DomainSpecificPlot(fun_name, args)

        setattr(parser.values, option.dest, out)

def view_file(filename, filter_names, options, view=None):
    if view is None:
        if options.show:
            offscreen = False

        else:
            offscreen = get_default(options.offscreen, True)
        view = Viewer(filename, watch=options.watch,
                      animate=options.anim_format is not None,
                      anim_format=options.anim_format,
                      ffmpeg_options=options.ffmpeg_options,
                      output_dir=options.output_dir,
                      offscreen=offscreen)

        if options.only_names is not None:
            options.only_names = options.only_names.split(',')

        view(show=options.show, is_3d=options.is_3d, view=options.view,
             roll=options.roll,
             fgcolor=options.fgcolor, bgcolor=options.bgcolor,
             layout=options.layout,
             scalar_mode=options.scalar_mode,
             vector_mode=options.vector_mode,
             rel_scaling=options.rel_scaling,
             clamping=options.clamping, ranges=options.ranges,
             is_scalar_bar=options.is_scalar_bar,
             is_wireframe=options.is_wireframe,
             opacity=options.opacity,
             subdomains_args=options.subdomains_args,
             rel_text_width=options.rel_text_width,
             fig_filename=options.filename, resolution=options.resolution,
             filter_names=filter_names, only_names=options.only_names,
             group_names=options.group_names,
             step=options.step, anti_aliasing=options.anti_aliasing,
             domain_specific=options.domain_specific)

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
    parser.add_option("", "--no-offscreen",
                      action="store_false", dest="offscreen",
                      default=None, help=help['no_offscreen'])
    parser.add_option("--3d",
                      action="store_true", dest="is_3d",
                      default=False, help=help['is_3d'])
    parser.add_option("--view", type='str',
                      metavar='angle,angle[,distance[,focal_point]]',
                      action="callback", dest="view",
                      callback=parse_view, help=help['view'])
    parser.add_option("--roll", type='float', metavar='angle',
                      action="store", dest="roll",
                      default=0.0, help=help['roll'])
    parser.add_option("--fgcolor", metavar='R,G,B',
                      action="store", dest="fgcolor",
                      default='0.0,0.0,0.0', help=help['fgcolor'])
    parser.add_option("--bgcolor", metavar='R,G,B',
                      action="store", dest="bgcolor",
                      default='1.0,1.0,1.0', help=help['bgcolor'])
    parser.add_option("--layout", metavar='layout',
                      action="store", dest="layout",
                      default='rowcol', help=help['layout'])
    parser.add_option("--scalar-mode", metavar='mode',
                      action="store", dest="scalar_mode",
                      default='iso_surface', help=help['scalar_mode'])
    parser.add_option("--vector-mode", metavar='mode',
                      action="store", dest="vector_mode",
                      default='arrows_norm', help=help['vector_mode'])
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
    parser.add_option("", "--wireframe",
                      action="store_true", dest="is_wireframe",
                      default=False, help=help['is_wireframe'])
    parser.add_option("--opacity", type='str', metavar='opacity',
                      action="callback", dest="opacity",
                      callback=parse_opacity, help=help['opacity'])
    parser.add_option("--rel-text-width", type='float', metavar='width',
                      action="store", dest="rel_text_width",
                      default=0.02, help=help['rel_text_width'])
    parser.add_option("-w", "--watch",
                      action="store_true", dest="watch",
                      default=False, help=help['watch'])
    parser.add_option("-o", "--output", metavar='filename',
                      action="store", dest="filename",
                      default=None, help=help['filename'])
    parser.add_option("--output-dir", metavar='directory',
                      action="store", dest="output_dir",
                      default=None, help=help['output_dir'])
    parser.add_option("-a", "--animation", metavar='<ffmpeg-supported format>',
                      action="store", dest="anim_format",
                      default=None, help=help['anim_format'])
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
    parser.add_option("--group-names", type='str', metavar='name1,...,nameN:...',
                      action="callback", dest="group_names",
                      callback=parse_group_names, help=help['group_names'])
    parser.add_option("--subdomains", type='str',
                      metavar='mat_id_name,threshold_limits,single_color',
                      action="callback", dest="subdomains_args",
                      callback=parse_subdomains, default=None,
                      help=help['subdomains'])
    parser.add_option("--step", type='int', metavar='step',
                      action="store", dest="step",
                      default=0, help=help['step'])
    parser.add_option("--anti-aliasing", type='int', metavar='value',
                      action="store", dest="anti_aliasing",
                      default=None, help=help['anti_aliasing'])
    parser.add_option("-d", "--domain-specific", type='str',
                      metavar="'var_name0,function_name0," \
                              "par0=val0,par1=val1,...:var_name1,...'",
                      action="callback", dest="domain_specific",
                      callback=parse_domain_specific, default=None,
                      help=help['domain_specific'])
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

    options.fgcolor = tuple([float(ii) for ii in
                             options.fgcolor.split(',')])
    assert_(len(options.fgcolor) == 3)

    options.bgcolor = tuple([float(ii) for ii in
                             options.bgcolor.split(',')])
    assert_(len(options.bgcolor) == 3)

    # Output dir / file names.
    if options.filename is None:
        options.filename = 'view.png'
        if options.output_dir is None:
            options.output_dir = '.'

    else:
        options.output_dir, options.filename = os.path.split(options.filename)

    # Data filtering,
    if not options.all:
        filter_names = ['node_groups', 'mat_id']
    else:
        filter_names = []

    if options.anim_format is not None:
        # Do not call show when saving an animation.
        options.show = False

    if options.list_ranges:
        all_ranges = {}
        for ii, filename in enumerate(filenames):
            output('%d: %s' % (ii, filename))

            file_source = create_file_source(filename)
            file_source.set_step(options.step)
            for key, val in get_data_ranges(file_source()).iteritems():
                all_ranges.setdefault(key, []).append(val[3:])
            
        if len(filenames) > 1:
            print 'summary of ranges:'
            for key, ranges in all_ranges.iteritems():
                aux = nm.array(ranges)
                print '  "%s": min(min): %s max(max): %s' % \
                      (key, aux[:,[0,2]].min(axis=0), aux[:,[1,3]].max(axis=0))

    else:
        if len(filenames) == 1:
            filenames = filenames[0]

        view = view_file(filenames, filter_names, options)        

    if options.anim_format is not None:
        view.encode_animation(options.filename, options.anim_format,
                              options.ffmpeg_options)

if __name__ == '__main__':
    main()

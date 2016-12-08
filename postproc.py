#!/usr/bin/env python
"""
This is a script for quick Mayavi-based visualizations of finite element
computations results.

Examples
--------
  The examples assume that run_tests.py has been run successfully and the
  resulting data files are present.

  - view data in output-tests/test_navier_stokes.vtk

    $ python postproc.py output-tests/test_navier_stokes.vtk
    $ python postproc.py output-tests/test_navier_stokes.vtk --3d

  - save a snapshot image and exit

    $ python postproc.py output-tests/test_poisson.vtk -o image.png -n

  - save a snapshot image without off-screen rendering and exit

    $ python postproc.py output-tests/test_poisson.vtk -o image.png \
                         -n --no-offscreen

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

    $ python postproc.py output-tests/test_hyperelastic_TL.*.vtk \
                         --ranges=u,0,0.02 -a mov --ffmpeg-options="-r 2 -sameq"
"""
from __future__ import print_function
from __future__ import absolute_import
from argparse import ArgumentParser, Action, RawDescriptionHelpFormatter
import os
import glob

import sfepy
from sfepy.base.base import assert_, get_default, output, nm
from sfepy.postprocess.viewer import (Viewer, get_data_ranges,
                                      create_file_source)
from sfepy.postprocess.domain_specific import DomainSpecificPlot
import six

helps = {
    'debug':
    'automatically start debugger when an exception is raised',
    'filename' :
    'view image file name [default: "view.png"]',
    'output_dir' :
    'output directory for saving view images; ignored when -o option is' \
    ' given, as the directory part of the filename is taken instead' \
    ' [default: "."]',
    'no_show' :
    'do not call mlab.show()',
    'no_offscreen' :
    'force no offscreen rendering for --no-show',
    'anim_format' :
    'if set to a ffmpeg-supported format (e.g. mov, avi, mpg), ffmpeg is' \
    ' installed and results of multiple time steps are given, an animation is' \
    ' created in the same directory as the view images',
    'ffmpeg_options' :
    'ffmpeg animation encoding options (enclose in "")' \
    '[default: "%(default)s"]',

    'step' :
    'set the time step. Negative indices are allowed, -1 means the last step.'
    ' The closest higher step is used if the desired one is not available.'
    ' Has precedence over --time. [default: the first step]',
    'time' :
    'set the time. The closest higher time is used if the desired one is not'
    ' available. [default: None]',
    'watch' :
    'watch the results file for changes (single file mode only)',
    'all' :
    'draw all data (normally, node_groups and mat_id are omitted)',
    'only_names' :
    'draw only named data',
    'list_ranges' :
    'do not plot, only list names and ranges of all data',
    'ranges' :
    'force data ranges [default: automatic from data]',

    'resolution' :
    'image resolution in NxN format [default: shorter axis: 600;'\
    ' depends on layout: for rowcol it is 800x600]',
    'layout' :
    'layout for multi-field plots, one of: rowcol, colrow, row, col, row#n,' \
    'col#n, where #n is the number of plots in the specified direction ' \
    '[default: %(default)s]',
    'is_3d' :
    '3d plot mode',
    'view' :
    'camera azimuth, elevation angles, and optionally also '
    'distance and focal point coordinates (without []) as in `mlab.view()` '
    '[default: if --3d is True: "45,45", else: "0,0"]',
    'roll' :
    'camera roll angle [default: %(default)s]',
    'parallel_projection' :
    'use parallel projection',
    'fgcolor' :
    'foreground color, that is the color of all text annotation labels'
    ' (axes, orientation axes, scalar bar labels) [default: %(default)s]',
    'bgcolor' :
    'background color [default: %(default)s]',
    'colormap' :
    'mayavi2 colormap name [default: %(default)s]',
    'anti_aliasing' :
    'value of anti-aliasing [default: mayavi2 default]',

    'is_scalar_bar' :
    'show scalar bar for each data',
    'is_wireframe' :
    'show wireframe of mesh surface for each data',
    'group_names' :
    'superimpose plots of data in each group',
    'subdomains' :
    'superimpose surfaces of subdomains over each data;' \
    ' example value: mat_id,0,None,True',
    'domain_specific' :
    'domain specific drawing functions and configurations',

    'scalar_mode' :
    'mode for plotting scalars with --3d, one of: cut_plane, iso_surface,'\
    ' both [default: %(default)s]',
    'vector_mode' :
    'mode for plotting vectors, one of: arrows, norm, arrows_norm, warp_norm'\
    ' [default: %(default)s]',
    'rel_scaling' :
    'relative scaling of glyphs (vector field visualization)' \
    ' [default: %(default)s]',
    'clamping' :
    'glyph clamping mode',
    'opacity' :
    'opacity in [0.0, 1.0]. Can be given either globally'
    ' as a single float, or per module, e.g.'
    ' "wireframe=0.1,scalar_cut_plane=0.5". Possible keywords are: wireframe,'
    ' scalar_cut_plane, vector_cut_plane, surface, iso_surface,'
    ' arrows_surface, glyphs. [default: 1.0]',
    'rel_text_width' :
    'relative text annotation width [default: %(default)s]',
}

class ParseView(Action):
    def __call__(self, parser, namespace, value, option_string=None):
        vals = value.split(',')
        assert_(len(vals) in [2, 3, 6])
        val = tuple(float(ii) for ii in vals)
        if len(vals) == 6:
            val = val[:3] + (list(val[3:]),)
        setattr(namespace, self.dest, val)

class ParseResolution(Action):
    def __call__(self, parser, namespace, value, option_string=None):
        if value is not None:
            print(value)
            setattr(namespace, self.dest,
                    tuple([int(r) for r in value.split('x')]))

class ParseRanges(Action):
    def __call__(self, parser, namespace, value, option_string=None):
        if value is not None:
            print(value)
            ranges = {}
            for rng in value.split(':'):
                aux = rng.split(',')
                ranges[aux[0]] = (float(aux[1]), float(aux[2]))
            setattr(namespace, self.dest, ranges)

class ParseOpacity(Action):
    def __call__(self, parser, namespace, value, option_string=None):
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

        setattr(namespace, self.dest, opacity)

class ParseGroupNames(Action):
    def __call__(self, parser, namespace, value, option_string=None):
        if value is not None:
            print(value)
            group_names = [tuple(group.split(','))
                           for group in value.split(':')]
            setattr(namespace, self.dest, group_names)

class ParseSubdomains(Action):
    def __call__(self, parser, namespace, value, option_string=None):
        if value is not None:
            print(value)
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
            setattr(namespace, self.dest, subdomains_args)

class ParseDomainSpecific(Action):
    def __call__(self, parser, namespace, value, option_string=None):
        if value is not None:
            print(value)
            out = {}
            confs = value.split(':')
            for conf in confs:
                aux = conf.split(',')
                var_name, fun_name = aux[:2]
                args = aux[2:]

                out[var_name] = DomainSpecificPlot(fun_name, args)
            setattr(namespace, self.dest, out)

def view_file(filename, filter_names, options, view=None):
    if view is None:
        if options.show:
            offscreen = False

        else:
            offscreen = get_default(options.offscreen, True)
        view = Viewer(filename, watch=options.watch,
                      ffmpeg_options=options.ffmpeg_options,
                      output_dir=options.output_dir,
                      offscreen=offscreen)

        if options.only_names is not None:
            options.only_names = options.only_names.split(',')

        view(show=options.show, is_3d=options.is_3d, view=options.view,
             roll=options.roll,
             parallel_projection=options.parallel_projection,
             fgcolor=options.fgcolor, bgcolor=options.bgcolor,
             colormap=options.colormap,
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
             step=options.step, time=options.time,
             anti_aliasing=options.anti_aliasing,
             domain_specific=options.domain_specific)

    else:
        view.set_source_filename(filename)
        view.save_image(options.filename)

    return view

def main():
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version',
                        version='%(prog)s ' + sfepy.__version__)
    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        default=False, help=helps['debug'])

    group = parser.add_argument_group('Output Options')
    group.add_argument('-o', '--output', metavar='filename',
                       action='store', dest='filename',
                       default=None, help=helps['filename'])
    group.add_argument('--output-dir', metavar='directory',
                       action='store', dest='output_dir',
                       default=None, help=helps['output_dir'])
    group.add_argument('-n', '--no-show',
                       action='store_false', dest='show',
                       default=True, help=helps['no_show'])
    group.add_argument('--no-offscreen',
                       action='store_false', dest='offscreen',
                       default=None, help=helps['no_offscreen'])
    group.add_argument('-a', '--animation',
                       metavar='<ffmpeg-supported format>', action='store',
                       dest='anim_format', default=None,
                       help=helps['anim_format'])
    group.add_argument('--ffmpeg-options', metavar='<ffmpeg options>',
                       action='store', dest='ffmpeg_options',
                       default='-r 10 -sameq',
                       help=helps['ffmpeg_options'])

    group = parser.add_argument_group('Data Options')
    group.add_argument('--step', type=int, metavar='step',
                       action='store', dest='step',
                       default=None, help=helps['step'])
    group.add_argument('--time', type=float, metavar='time',
                       action='store', dest='time',
                       default=None, help=helps['time'])
    group.add_argument('-w', '--watch',
                       action='store_true', dest='watch',
                       default=False, help=helps['watch'])
    group.add_argument('--all',
                       action='store_true', dest='all',
                       default=False, help=helps['all'])
    group.add_argument('--only-names', metavar='list of names',
                       action='store', dest='only_names',
                       default=None, help=helps['only_names'])
    group.add_argument('-l', '--list-ranges',
                       action='store_true', dest='list_ranges',
                       default=False, help=helps['list_ranges'])
    group.add_argument('--ranges', type=str,
                       metavar='name1,min1,max1:name2,min2,max2:...',
                       action=ParseRanges, dest='ranges',
                       help=helps['ranges'])

    group = parser.add_argument_group('View Options')
    group.add_argument('-r', '--resolution', type=str, metavar='resolution',
                       action=ParseResolution, dest='resolution',
                       help=helps['resolution'])
    group.add_argument('--layout', metavar='layout',
                       action='store', dest='layout',
                       default='rowcol', help=helps['layout'])
    group.add_argument('--3d',
                       action='store_true', dest='is_3d',
                       default=False, help=helps['is_3d'])
    group.add_argument('--view', type=str,
                       metavar='angle,angle[,distance[,focal_point]]',
                       action=ParseView, dest='view',
                       help=helps['view'])
    group.add_argument('--roll', type=float, metavar='angle',
                       action='store', dest='roll',
                       default=0.0, help=helps['roll'])
    group.add_argument('--parallel-projection',
                       action='store_true', dest='parallel_projection',
                       default=False, help=helps['parallel_projection'])
    group.add_argument('--fgcolor', metavar='R,G,B',
                       action='store', dest='fgcolor',
                       default='0.0,0.0,0.0', help=helps['fgcolor'])
    group.add_argument('--bgcolor', metavar='R,G,B',
                       action='store', dest='bgcolor',
                       default='1.0,1.0,1.0', help=helps['bgcolor'])
    group.add_argument('--colormap', metavar='colormap',
                       action='store', dest='colormap',
                       default='blue-red', help=helps['colormap'])
    group.add_argument('--anti-aliasing', type=int, metavar='value',
                       action='store', dest='anti_aliasing',
                       default=None, help=helps['anti_aliasing'])

    group = parser.add_argument_group('Custom Plots Options')
    group.add_argument('-b', '--scalar-bar',
                       action='store_true', dest='is_scalar_bar',
                       default=False, help=helps['is_scalar_bar'])
    group.add_argument('--wireframe',
                       action='store_true', dest='is_wireframe',
                       default=False, help=helps['is_wireframe'])
    group.add_argument('--group-names', type=str,
                       metavar='name1,...,nameN:...', action=ParseGroupNames,
                       dest='group_names', help=helps['group_names'])
    group.add_argument('--subdomains', type=str,
                       metavar='mat_id_name,threshold_limits,single_color',
                       action=ParseSubdomains, dest='subdomains_args',
                       default=None,
                       help=helps['subdomains'])
    group.add_argument('-d', '--domain-specific', type=str,
                       metavar='"var_name0,function_name0,' \
                       'par0=val0,par1=val1,...:var_name1,..."',
                       action=ParseDomainSpecific, dest='domain_specific',
                       default=None,
                       help=helps['domain_specific'])

    group = parser.add_argument_group('Mayavi Options')
    group.add_argument('--scalar-mode', metavar='mode',
                       action='store', dest='scalar_mode',
                       default='iso_surface', help=helps['scalar_mode'])
    group.add_argument('--vector-mode', metavar='mode',
                       action='store', dest='vector_mode',
                       default='arrows_norm', help=helps['vector_mode'])
    group.add_argument('-s', '--scale-glyphs', type=float, metavar='scale',
                       action='store', dest='rel_scaling',
                       default=0.05, help=helps['rel_scaling'])
    group.add_argument('--clamping',
                       action='store_true', dest='clamping',
                       default=False, help=helps['clamping'])
    group.add_argument('--opacity', type=str, metavar='opacity',
                       action=ParseOpacity, dest='opacity',
                       help=helps['opacity'])
    group.add_argument('--rel-text-width', type=float, metavar='width',
                       action='store', dest='rel_text_width',
                       default=0.02, help=helps['rel_text_width'])

    parser.add_argument('filenames', nargs='+')
    options = parser.parse_args()

    if options.debug:
        from sfepy.base.base import debug_on_error; debug_on_error()

    filenames = options.filenames

    options.fgcolor = tuple([float(ii) for ii in
                             options.fgcolor.split(',')])
    assert_(len(options.fgcolor) == 3)

    options.bgcolor = tuple([float(ii) for ii in
                             options.bgcolor.split(',')])
    assert_(len(options.bgcolor) == 3)

    can_save = not options.show

    # Output dir / file names.
    if options.filename is None:
        can_save = False
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
            if (options.step is None) and (options.time is None):
                steps, _ = file_source.get_ts_info()

            else:
                if options.step is not None:
                    step, _ = file_source.get_step_time(step=options.step)

                else:
                    step, _ = file_source.get_step_time(time=options.time)

                steps = [step]

            if not len(steps):
                steps = [0]

            for iis, step in enumerate(steps):
                output('%d: step %d' %(iis, step))
                file_source.get_step_time(step=step)
                source = file_source.create_source()
                ranges = get_data_ranges(source, return_only=True)
                for key, val in six.iteritems(ranges):
                    all_ranges.setdefault(key, []).append(val[3:])

        if (len(filenames) > 1) or (len(steps) > 1):
            output('union of ranges:')

        else:
            output('ranges:')

        for key, ranges in six.iteritems(all_ranges):
            aux = nm.array(ranges)
            mins = aux[:, [0, 2]].min(axis=0)
            maxs = aux[:, [1, 3]].max(axis=0)
            output('  items: %s,%e,%e' % (key, mins[0], maxs[0]))
            output('  norms: %s,%e,%e' % (key, mins[1], maxs[1]))

    else:
        if len(filenames) == 1:
            filenames = filenames[0]

        view = view_file(filenames, filter_names, options)
        if can_save:
            view.save_image(options.filename)

    if options.anim_format is not None:
        view.save_animation(options.filename)
        view.encode_animation(options.filename, options.anim_format,
                              options.ffmpeg_options)

if __name__ == '__main__':
    main()

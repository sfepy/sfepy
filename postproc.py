#!/usr/bin/env python
from optparse import OptionParser

import sfepy
from sfepy.base.base import pause
from sfepy.postprocess import Viewer

usage = """%prog [options] filename"""

help = {
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
    'rel_text_width' :
    'relative text annotation width [default: %default]',
    'filename' :
    'view image file name' \
    ' [default: %default]',
    'all' :
    'draw all data (normally, node_groups and mat_id are omitted)',
    'list_names' :
    'do not plot, only list all dataset names',
    'only_names' :
    'draw only named data',
}

def parse_view(option, opt, value, parser):
    val = tuple(float(ii) for ii in value.split(','))
    setattr(parser.values, option.dest, val)


def main():
    parser = OptionParser(usage=usage, version="%prog " + sfepy.__version__)
    parser.add_option("--no-show",
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
    parser.add_option("--rel-text-width", type='float', metavar='width',
                      action="store", dest="rel_text_width",
                      default=0.02, help=help['rel_text_width'])
    parser.add_option("-o", metavar='filename',
                      action="store", dest="filename",
                      default='view.png', help=help['filename'])
    parser.add_option("-a", "--all",
                      action="store_true", dest="all",
                      default=False, help=help['all'])
    parser.add_option("--list-names",
                      action="store_true", dest="list_names",
                      default=False, help=help['list_names'])
    parser.add_option("--only-names", metavar='list of names',
                      action="store", dest="only_names",
                      default=None, help=help['only_names'])
    options, args = parser.parse_args()

    if (len(args) == 1):
        filename = args[0]
    else:
        parser.print_help(),
        return

    if not options.all:
        filter_names = ['node_groups', 'mat_id']
    else:
        filter_names = []

    view = Viewer(filename, offscreen=not options.show)

    if options.list_names:
        for line in view.get_data_names():
            print line

    else:
        if options.only_names is not None:
            options.only_names = options.only_names.split(',')

        view(show=options.show, is_3d=options.is_3d, view=options.view,
             roll=options.roll, layout=options.layout,
             scalar_mode=options.scalar_mode, rel_scaling=options.rel_scaling,
             clamping=options.clamping, rel_text_width=options.rel_text_width,
             fig_filename=options.filename, filter_names=filter_names,
             only_names=options.only_names)

if __name__ == '__main__':
    main()

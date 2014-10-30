#!/usr/bin/env python
"""
Plot quadrature points for the given geometry and integration order.
"""
import sys
sys.path.append('.')
from optparse import OptionParser

import sfepy.postprocess.plot_quadrature as pq

usage = '%prog [options]\n' + __doc__.rstrip()

helps = {
    'geometry' :
    'reference element geometry, one of "2_3", "2_4", "3_4", "3_8"'
    ' [default: %default]',
    'order' :
    'quadrature order [default: %default]',
    'min_radius' :
    'min. radius of points corresponding to the min. weight'
    ' [default:  %default]',
    'max_radius' :
    'max. radius of points corresponding to the max. weight'
    ' [default:  %default]',
    'show_colorbar' :
    'show colorbar for quadrature weights'
}

def main():
    parser = OptionParser(usage=usage, version='%prog')
    parser.add_option('-g', '--geometry', metavar='name',
                      action='store', dest='geometry',
                      default='2_4', help=helps['geometry'])
    parser.add_option('-n', '--order', metavar='order', type=int,
                      action='store', dest='order',
                      default=2, help=helps['order'])
    parser.add_option('-r', '--min-radius', metavar='float', type=float,
                      action='store', dest='min_radius',
                      default=10, help=helps['min_radius'])
    parser.add_option('-R', '--max-radius', metavar='float', type=float,
                      action='store', dest='max_radius',
                      default=50, help=helps['max_radius'])
    parser.add_option('-c', '--show-colorbar',
                      action='store_true', dest='show_colorbar',
                      default=False, help=helps['show_colorbar'])
    options, args = parser.parse_args()

    if len(args) != 0:
        parser.print_help(),
        return

    pq.plot_quadrature(None, options.geometry, options.order,
                       options.min_radius, options.max_radius,
                       options.show_colorbar)
    pq.plt.show()

if __name__ == '__main__':
    main()

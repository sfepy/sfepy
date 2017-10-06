#!/usr/bin/env python
"""
Plot quadrature points for the given geometry and integration order.
"""
from __future__ import absolute_import, print_function
import sys
sys.path.append('.')
from argparse import ArgumentParser

import sfepy.postprocess.plot_quadrature as pq

helps = {
    'geometry' :
    'reference element geometry, one of "2_3", "2_4", "3_4", "3_8"'
    ' [default: %(default)s]',
    'order' :
    'quadrature order [default: %(default)s]',
    'boundary' :
    'plot boundary quadrature points',
    'min_radius' :
    'min. radius of points corresponding to the min. weight'
    ' [default:  %(default)s]',
    'max_radius' :
    'max. radius of points corresponding to the max. weight'
    ' [default:  %(default)s]',
    'show_colorbar' :
    'show colorbar for quadrature weights',
    'show_labels' :
    'label quadrature points',
    'print_qp' :
    'print quadrature points and weights',
}

def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version='%(prog)s')
    parser.add_argument('-g', '--geometry', metavar='name',
                        action='store', dest='geometry',
                        default='2_4', help=helps['geometry'])
    parser.add_argument('-n', '--order', metavar='order', type=int,
                        action='store', dest='order',
                        default=2, help=helps['order'])
    parser.add_argument('-b', '--boundary',
                        action='store_true', dest='boundary',
                        default=False, help=helps['boundary'])
    parser.add_argument('-r', '--min-radius', metavar='float', type=float,
                        action='store', dest='min_radius',
                        default=10, help=helps['min_radius'])
    parser.add_argument('-R', '--max-radius', metavar='float', type=float,
                        action='store', dest='max_radius',
                        default=50, help=helps['max_radius'])
    parser.add_argument('-c', '--show-colorbar',
                        action='store_true', dest='show_colorbar',
                        default=False, help=helps['show_colorbar'])
    parser.add_argument('-l', '---show-labels',
                        action='store_true', dest='show_labels',
                        default=False, help=helps['show_labels'])
    parser.add_argument('-p', '--print-qp',
                        action='store_true', dest='print_qp',
                        default=False, help=helps['print_qp'])
    options = parser.parse_args()

    aux = pq.plot_quadrature(None, options.geometry, options.order,
                             boundary=options.boundary,
                             min_radius=options.min_radius,
                             max_radius=options.max_radius,
                             show_colorbar=options.show_colorbar,
                             show_labels=options.show_labels)

    if options.print_qp:
        ax, coors, weights = aux
        for ic, coor in enumerate(coors):
            print(ic, coor, weights[ic])

    pq.plt.show()

if __name__ == '__main__':
    main()

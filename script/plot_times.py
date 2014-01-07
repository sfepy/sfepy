#!/usr/bin/env python
"""
Plot time steps, times of time steps and time deltas in a HDF5 results file.
"""
from optparse import OptionParser

import numpy as nm
import matplotlib.pyplot as plt

from sfepy.postprocess.time_history import extract_times

usage = '%prog [options] filename\n' + __doc__.rstrip()

helps = {
    'logarithmic' :
    'plot time steps in logarithmic scale',
}

def main():
    parser = OptionParser(usage=usage, version='%prog')
    parser.add_option('-l', '--logarithmic',
                      action='store_true', dest='logarithmic',
                      default=False, help=helps['logarithmic'])
    options, args = parser.parse_args()

    if (len(args) == 1):
        filename = args[0]
    else:
        parser.print_help()
        return

    plt.rcParams['lines.linewidth'] = 3
    plt.rcParams['lines.markersize'] = 9
    fontsize = 16

    steps, times, nts, dts = extract_times(filename)
    dts[-1] = nm.nan

    ax = plt.subplot(211)
    if options.logarithmic:
        l1, = ax.semilogy(steps, dts, 'b')
    else:
        l1, = ax.plot(steps, dts, 'b')
    ax.set_xlabel('step', fontsize=fontsize)
    ax.set_ylabel(r'$\Delta t$', fontsize=fontsize)
    ax.grid(True)

    ax = ax.twinx()
    l2, = ax.plot(steps, times, 'g')
    ax.set_ylabel(r'$t$', fontsize=fontsize)
    ax.legend([l1, l2], [r'$\Delta t$', r'$t$'], loc=0)

    ax = plt.subplot(212)
    if options.logarithmic:
        ax.semilogy(times, dts, 'b+')
    else:
        ax.plot(times, dts, 'b+')
    ax.set_xlabel(r'$t$', fontsize=fontsize)
    ax.set_ylabel(r'$\Delta t$', fontsize=fontsize)
    ax.grid(True)

    plt.show()

if __name__ == '__main__':
    main()

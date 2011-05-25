#!/usr/bin/env python
"""
Simple script for testing various SfePy functionality, examples not
covered by tests, and running the tests.

The script just runs the commands specified in its main() using the
`subprocess` module, captures the output and compares one or more key
words to the expected ones.

The output of failed commands is saved to 'test_install.log' file.
"""
import time
from optparse import OptionParser
import shlex
import subprocess

def check_output(cmd):
    """
    Run the specified command and capture its outputs.

    Returns
    -------
    out : tuple
        The (stdout, stderr) output tuple.
    """
    print cmd
    args = shlex.split(cmd)

    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = p.communicate()

    return out

def report(out, name, line, item, value, eps=None):
    """
    Check that `item` at `line` of the output string `out` is equal
    to `value`. If not, print the output.
    """
    try:
        status = out.split('\n')[line].split()

    except IndexError:
        print '  not enough output from command!'
        ok = False

    else:
        try:
            print '  comparing:', status[item], value

            if eps is None:
                ok = (status[item] == value)

            else:
                try:
                    ok = abs(float(status[item]) - float(value)) < eps

                except:
                    ok = False

        except IndexError:
            ok = False

    if ok:
        print '  %s:' % name, ok

    else:
        print '! %s:' % name, ok

        fd = open('test_install.log', 'a')
        fd.write('*' * 55)
        fd.write(out)
        fd.write('*' * 55)

usage = '%prog' + '\n' + __doc__

def main():
    parser = OptionParser(usage=usage, version='%prog')
    options, args = parser.parse_args()
    if len(args) > 0:
        parser.print_help()
        return

    fd = open('test_install.log', 'w')
    fd.close()

    t0 = time.time()

    out, err = check_output('./simple.py examples/diffusion/poisson.py')
    report(out, '...', -2, 5, '1.173819e-16', eps=1e-15)

    out, err = check_output('./simple.py examples/navier_stokes/stokes.py')
    report(out, '...', -2, 5, '1.210678e-13', eps=1e-11)

    out, err = check_output('./simple.py examples/diffusion/poisson_parametric_study.py')
    report(out, '...', -2, 5, '1.606408e-14', eps=1e-13)

    out, err = check_output('./findSurf.py meshes/quantum/cube.node -')
    report(out, '...', -2, 1, '64247')

    out, err = check_output('./eigen.py examples/phononic/band_gaps.py')
    report(out, '...', -4, 2, '232.40156299')
    report(out, '...', -3, 1, '132604.79235405]')

    out, err = check_output('./eigen.py examples/phononic/band_gaps.py -d')
    report(out, '...', -5, 4, '0.209329,')

    out, err = check_output('./schroedinger.py --2d --mesh')
    report(out, '...', -2, -1, 'tmp/mesh.vtk')

    out, err = check_output('./schroedinger.py --2d --hydrogen')
    report(out, '...', -4, -2, '-0.01984415', eps=1e-4)

    out, err = check_output('python examples/standalone/homogenized_elasticity/rs_correctors.py -n')
    report(out, '...', -2, -1, '1.644e-01]]')

    out, err = check_output('python examples/standalone/elastic_materials/compare_elastic_materials.py -n')
    report(out, '...', -2, 5, '1.068759e-14', eps=1e-13)

    out, err = check_output('python examples/standalone/interactive/linear_elasticity.py')
    report(out, '...', -8, 0, '1.62128841139e-14', eps=1e-13)

    t1 = time.time()

    out, err = check_output('./runTests.py')
    report(out, 'tests', -2, 7, '0')

    t2 = time.time()

    fd = open('test_install_times.log', 'a+')
    fd.write('%s: examples: %.2f [s], tests: %.2f [s]\n'
             % (time.ctime(t0), t1 - t0, t2 - t1))
    fd.close()

if __name__ == '__main__':
    main()

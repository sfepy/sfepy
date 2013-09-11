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

def report(out, name, line, item, value, eps=None, return_item=False):
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

    if return_item:
        return ok, status[item]

    else:
        return ok

usage = '%prog' + '\n' + __doc__

def main():
    parser = OptionParser(usage=usage, version='%prog')
    options, args = parser.parse_args()
    if len(args) > 0:
        parser.print_help()
        return

    fd = open('test_install.log', 'w')
    fd.close()

    eok = 0

    t0 = time.time()

    out, err = check_output('python ./simple.py examples/diffusion/poisson.py')
    eok += report(out, '...', -2, 5, '1.173819e-16', eps=1e-15)

    out, err = check_output("""python ./simple.py -c "ebc_2 : {'name' : 't2', 'region' : 'Gamma_Right', 'dofs' : {'t.0' : -5.0}}" examples/diffusion/poisson.py""")
    eok += report(out, '...', -2, 5, '2.308051e-16', eps=1e-15)

    out, err = check_output('python ./simple.py examples/navier_stokes/stokes.py')
    eok += report(out, '...', -2, 5, '1.210678e-13', eps=1e-11)

    out, err = check_output('python ./simple.py examples/diffusion/poisson_parametric_study.py')
    eok += report(out, '...', -2, 5, '1.606408e-14', eps=1e-13)

    out, err = check_output('python ./simple.py examples/linear_elasticity/its2D_3.py')
    eok += report(out, '...', -23, 5, '3.964886e-12', eps=1e-11)
    eok += report(out, '...', -3, 4, '2.58660e+01', eps=1e-5)

    out, err = check_output('python ./simple.py examples/linear_elasticity/linear_elastic_probes.py')
    eok += report(out, '...', -11, 5, '4.638192e-18', eps=1e-15)

    out, err = check_output('python ./findSurf.py meshes/quantum/cube.node -')
    eok += report(out, '...', -2, 0, '64247')

    out, err = check_output('python ./phonon.py examples/phononic/band_gaps.py')
    eok += report(out, '...', -6, 2, '208.54511594')
    eok += report(out, '...', -5, 1, '116309.22337295]')

    out, err = check_output('python ./phonon.py examples/phononic/band_gaps.py --phase-velocity')
    eok += report(out, '...', -2, 3, '4.1894123')
    eok += report(out, '...', -2, 4, '2.62055608]')

    out, err = check_output('python ./phonon.py examples/phononic/band_gaps.py -d')
    eok += report(out, '...', -6, 1, '[0,')

    out, err = check_output('python ./phonon.py examples/phononic/band_gaps_rigid.py')
    eok += report(out, '...', -6, 2, '4.58709531e+01')
    eok += report(out, '...', -5, 1, '1.13929200e+05]')

    out, err = check_output('python ./schroedinger.py --2d --create-mesh')
    eok += report(out, '...', -2, -1, 'tmp/mesh.vtk')

    out, err = check_output('python ./schroedinger.py --2d --hydrogen')
    eok += report(out, '...', -4, -2, '-0.01984415', eps=1e-4)

    out, err = check_output('python ./homogen.py examples/homogenization/perfusion_micro.py')
    eok += report(out, '...', -7, -1, 'EpA...')

    out, err = check_output('python examples/standalone/homogenized_elasticity/rs_correctors.py -n')
    eok += report(out, '...', -2, -1, '1.644e-01]]')

    out, err = check_output('python examples/standalone/elastic_materials/compare_elastic_materials.py -n')
    eok += report(out, '...', -2, 5, '1.068759e-14', eps=1e-13)

    out, err = check_output('python examples/standalone/interactive/linear_elasticity.py')
    eok += report(out, '...', -8, 0, '1.62128841139e-14', eps=1e-13)

    out, err = check_output('python examples/standalone/thermal_electric/thermal_electric.py')
    eok += report(out, '...', -3, 5, '2.612933e-14', eps=1e-13)

    t1 = time.time()

    out, err = check_output('python ./runTests.py')
    tok, failed = report(out, 'tests', -2, 7, '0', return_item=True)
    tok = {True : 'ok', False : 'fail'}[tok]

    t2 = time.time()

    fd = open('test_install_times.log', 'a+')
    fd.write('%s: examples: %.2f [s] (%d), tests: %.2f [s] (%s: %s)\n'
             % (time.ctime(t0), t1 - t0, eok, t2 - t1, tok, failed))
    fd.close()

if __name__ == '__main__':
    main()

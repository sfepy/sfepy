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

def report2(out, name, items, return_item=False):
    """
    Check that `items` are in the output string `out`.
    If not, print the output.
    """
    ok = True
    for s in items:
        print '  checking:', s
        if s not in out:
            ok = False
            break

    if ok:
        print '  %s:' % name, ok

    else:
        print '! %s:' % name, ok

        fd = open('test_install.log', 'a')
        fd.write('*' * 55)
        fd.write(out)
        fd.write('*' * 55)

    if return_item:
        return ok, s

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

    out, err = check_output('python ./script/blockgen.py')
    eok += report(out, '...', -2, 1, '...done')

    out, err = check_output('python ./script/cylindergen.py')
    eok += report(out, '...', -2, 1, '...done')

    out, err = check_output('python ./script/convert_mesh.py meshes/3d/cylinder.vtk out.mesh')
    eok += report(out, '...', -2, 1, '...done')

    out, err = check_output('python ./script/tile_periodic_mesh.py -r 2,2 meshes/elements/2_4_2.mesh out-per.mesh')
    eok += report(out, '...', -2, 1, 'done.')

    out, err = check_output('python ./script/extract_surface.py meshes/various_formats/octahedron.node -')
    eok += report(out, '...', -2, 0, '1185')

    out, err = check_output('python ./simple.py examples/diffusion/poisson.py')
    eok += report(out, '...', -2, 5, '1.173819e-16', eps=1e-15)

    out, err = check_output("""python ./simple.py -c "ebc_2 : {'name' : 't2', 'region' : 'Gamma_Right', 'dofs' : {'t.0' : -5.0}}" examples/diffusion/poisson.py""")
    eok += report(out, '...', -2, 5, '2.308051e-16', eps=1e-15)

    out, err = check_output('python ./simple.py examples/diffusion/poisson_iga.py')
    eok += report(out, '...', -2, 5, '3.373487e-15', eps=1e-14)

    out, err = check_output('python ./simple.py examples/navier_stokes/stokes.py')
    eok += report(out, '...', -2, 5, '1.210678e-13', eps=1e-11)

    out, err = check_output('python ./simple.py examples/diffusion/poisson_parametric_study.py')
    eok += report(out, '...', -2, 5, '1.606408e-14', eps=1e-13)

    out, err = check_output('python ./simple.py examples/linear_elasticity/its2D_3.py')
    eok += report(out, '...', -23, 5, '3.964886e-12', eps=1e-11)
    eok += report(out, '...', -3, 4, '2.58660e+01', eps=1e-5)

    out, err = check_output('python ./simple.py examples/linear_elasticity/linear_elastic.py --format h5')
    eok += report(out, '...', -2, 5, '4.638192e-18', eps=1e-15)

    out, err = check_output('python ./extractor_sfepy.py -d cylinder.h5')
    eok += report(out, '...', -2, 1, '...done')

    out, err = check_output('python ./postproc.py -n --no-offscreen -o cylinder.png cylinder.h5')
    eok += report(out, '...', -3, 2, 'cylinder.png...')

    out, err = check_output('python ./phonon.py examples/phononic/band_gaps.py')
    eok += report(out, '...', -7, 2, '208.54511594')
    eok += report(out, '...', -6, 1, '116309.22337295]')

    out, err = check_output('python ./phonon.py examples/phononic/band_gaps.py --phase-velocity')
    eok += report(out, '...', -2, 3, '4.1894123')
    eok += report(out, '...', -2, 4, '2.62055608]')

    out, err = check_output('python ./phonon.py examples/phononic/band_gaps.py -d')
    eok += report(out, '...', -6, 1, '[0,')

    out, err = check_output('python ./phonon.py examples/phononic/band_gaps_rigid.py')
    eok += report(out, '...', -7, 2, '4.58709531e+01')
    eok += report(out, '...', -6, 1, '1.13929200e+05]')

    out, err = check_output('python ./schroedinger.py --hydrogen')
    eok += report(out, '...', -4, -2, '-0.01913506', eps=1e-4)

    out, err = check_output('python ./homogen.py examples/homogenization/perfusion_micro.py')
    eok += report2(out, '...', ['computing EpA', 'computing PA_3',
                                'computing GA', 'computing EmA',
                                'computing KA'])

    out, err = check_output('python examples/homogenization/rs_correctors.py -n')
    eok += report(out, '...', -2, -1, '1.644e-01]]')

    out, err = check_output('python examples/large_deformation/compare_elastic_materials.py -n')
    eok += report(out, '...', -2, 5, '1.068759e-14', eps=1e-13)

    out, err = check_output('python examples/linear_elasticity/linear_elastic_interactive.py')
    eok += report(out, '...', -8, 0, '1.62128841139e-14', eps=1e-13)

    out, err = check_output('python examples/linear_elasticity/modal_analysis.py')
    eok += report(out, '...', -12, 5, '12142.11470773', eps=1e-13)

    out, err = check_output('python examples/multi_physics/thermal_electric.py')
    eok += report(out, '...', -3, 5, '2.612933e-14', eps=1e-13)

    out, err = check_output('mpiexec -n 2 python examples/diffusion/poisson_parallel_interactive.py output-parallel -2 --silent -ksp_monitor')
    eok += report(out, '...', -2, 4, '8.021313824020e-07', eps=1e-6)

    out, err = check_output('mpiexec -n 2 python examples/multi_physics/biot_parallel_interactive.py output-parallel -2 --silent -ksp_monitor')
    eok += report(out, '...', -2, 4, '3.787214380277e-09', eps=1e-8)

    t1 = time.time()

    out, err = check_output('python ./run_tests.py')
    tok, failed = report(out, 'tests', -2, 7, '0', return_item=True)
    tok = {True : 'ok', False : 'fail'}[tok]

    t2 = time.time()

    fd = open('test_install_times.log', 'a+')
    fd.write('%s: examples: %.2f [s] (%d), tests: %.2f [s] (%s: %s)\n'
             % (time.ctime(t0), t1 - t0, eok, t2 - t1, tok, failed))
    fd.close()

if __name__ == '__main__':
    main()

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
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import shlex
import subprocess
import logging
import re

DEBUG_FMT = '*' * 55 + '\n%s\n' + '*' * 55

def _get_logger(filename='test_install.log'):
    """
    Convenience function to set-up output and logging.
    """
    logger = logging.getLogger('test_install.py')
    logger.setLevel(logging.DEBUG)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)

    file_handler = logging.FileHandler(filename)
    file_handler.setLevel(logging.DEBUG)

    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    return logger

logger = _get_logger()

def check_output(cmd):
    """
    Run the specified command and capture its outputs.

    Returns
    -------
    out : tuple
        The (stdout, stderr) output tuple.
    """
    logger.info(cmd)
    args = shlex.split(cmd)

    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = [ii.decode() for ii in p.communicate()]

    return out

def report(out, name, line, item, value, eps=None, return_item=False,
           match_numbers=False):
    """
    Check that `item` at `line` of the output string `out` is equal
    to `value`. If not, print the output.
    """
    try:
        if match_numbers:
            status = out.split('\n')[line]

        else:
            status = out.split('\n')[line].split()

    except IndexError:
        logger.error('  not enough output from command!')
        ok = False

    else:
        try:
            if match_numbers:
                pat = r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?[jJ]?)'
                matches = re.findall(pat, status)
                status_item = matches[item]

            else:
                status_item = status[item]

            logger.info('  comparing: %s %s', status_item, value)

            if eps is None:
                ok = (status_item == value)

            else:
                try:
                    ok = abs(float(status_item) - float(value)) < eps

                except:
                    ok = False

        except IndexError:
            ok = False

    logger.info('  %s: %s', name, ok)

    if not ok:
        logger.debug(DEBUG_FMT, out)

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
        logger.info('  checking: %s', s)
        if s not in out:
            ok = False
            break

    logger.info('  %s: %s', name, ok)

    if not ok:
        logger.debug(DEBUG_FMT, out)

    if return_item:
        return ok, s

    else:
        return ok

def report_tests(out, return_item=False):
    """
    Check that all tests in the output string `out` passed.
    If not, print the output.
    """
    from pyparsing import (Word, Combine, Suppress, Optional, OneOrMore,
                           delimitedList, nums, Literal)
    from functools import partial
    integer = Word(nums).setName('integer')
    real = Combine(Word(nums) + '.' + Optional(Word(nums))).setName('real')

    equals = Suppress(OneOrMore('='))
    _stats = {}
    def add_stat(s, loc, toks, key=None):
        if key is None:
            key = toks[1]
        _stats[key] = toks[0]
        return toks

    word = ((integer + 'failed') |
            (integer + 'passed') |
            (integer + 'deselected') |
            (integer + 'warnings')).setParseAction(add_stat)

    line = (equals +
            Optional(delimitedList(word)) +
            'in' +
            (real + (Literal('s') | 'seconds'))
            .setParseAction(partial(add_stat, key='seconds')) +
            equals)

    line.searchString(out)

    keys = ['failed', 'passed', 'deselected', 'warnings', 'seconds']
    stats = {key : _stats.get(key, '0') for key in keys}

    ok = stats['failed'] == '0'

    logger.info(
        ('  {failed} failed, {passed} passed, {deselected} deselected,'
         ' {warnings} warnings in {seconds} seconds').format(**stats)
    )

    if not ok:
        logger.debug(DEBUG_FMT, out)

    if return_item:
        return ok, stats['failed']

    else:
        return ok

def main():
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s')
    parser.parse_args()

    fd = open('test_install.log', 'w')
    fd.close()

    eok = 0

    t0 = time.time()

    out, err = check_output('python3 sfepy/scripts/blockgen.py')
    eok += report(out, '...', -2, 1, '...done')

    out, err = check_output('python3 sfepy/scripts/cylindergen.py')
    eok += report(out, '...', -2, 1, '...done')

    out, err = check_output('python3 sfepy/scripts/convert_mesh.py meshes/3d/cylinder.vtk out.mesh')
    eok += report(out, '...', -2, 1, '...done')

    out, err = check_output('python3 sfepy/scripts/convert_mesh.py --tile 2,2 meshes/elements/2_4_2.mesh out-per.mesh')
    eok += report(out, '...', -2, 1, '...done')

    out, err = check_output('python3 sfepy/scripts/convert_mesh.py --extract-surface --print-surface=- meshes/various_formats/octahedron.node surf_octahedron.mesh')
    eok += report(out, '...', -4, 0, '1185')

    out, err = check_output("""python3 sfepy/scripts/simple.py -c "ebc_2 : {'name' : 't2', 'region' : 'Gamma_Right', 'dofs' : {'t.0' : -5.0}}" sfepy/examples/diffusion/poisson.py""")
    eok += report(out, '...', -8, 5, '2.308051e-16', eps=1e-15)

    out, err = check_output('python3 sfepy/scripts/simple.py sfepy/examples/diffusion/poisson_parametric_study.py')
    eok += report(out, '...', -8, 5, '1.606408e-14', eps=1e-13)

    out, err = check_output('python3 sfepy/scripts/simple.py sfepy/examples/linear_elasticity/its2D_3.py')
    eok += report(out, '...', -29, 5, '3.964886e-12', eps=1e-11)
    eok += report(out, '...', -9, 4, '2.58660e+01', eps=1e-5)

    out, err = check_output('python3 sfepy/scripts/simple.py sfepy/examples/linear_elasticity/linear_elastic.py --format h5')
    eok += report(out, '...', -8, 5, '4.638192e-18', eps=1e-15)

    out, err = check_output('python3 sfepy/scripts/extractor.py -d cylinder.h5')
    eok += report(out, '...', -2, 1, '...done')

    out, err = check_output('python3 sfepy/scripts/resview.py --off-screen -o cylinder.png cylinder.h5')
    eok += report(out, '...', -2, 1, 'cylinder.png')

    out, err = check_output('python3 sfepy/scripts/simple.py sfepy/examples/phononic/band_gaps.py')
    eok += report(out, '...', -9, 0, '2.08545116e+08', match_numbers=True)
    eok += report(out, '...', -8, 1, '1.16309223e+11', match_numbers=True)

    out, err = check_output('python3 sfepy/scripts/simple.py sfepy/examples/phononic/band_gaps.py --phonon-phase-velocity')
    eok += report(out, '...', -2, 0, '4189.41229592', match_numbers=True)
    eok += report(out, '...', -2, 1, '2620.55608256', match_numbers=True)

    out, err = check_output('python3 sfepy/scripts/simple.py sfepy/examples/phononic/band_gaps.py --phonon-dispersion')
    eok += report(out, '...', -6, 1, '[0,')

    out, err = check_output('python3 sfepy/scripts/simple.py sfepy/examples/phononic/band_gaps_rigid.py')
    eok += report(out, '...', -9, 0, '4.58709531e+07', match_numbers=True)
    eok += report(out, '...', -8, 1, '1.13929200e+11', match_numbers=True)

    out, err = check_output('python3 sfepy/scripts/simple.py sfepy/examples/quantum/hydrogen.py')
    eok += report(out, '...', -2, -2, '-0.01913506', eps=1e-4)

    out, err = check_output('python3 sfepy/scripts/simple.py sfepy/examples/homogenization/perfusion_micro.py')
    eok += report2(out, '...', ['computing EpA', 'computing PA_3',
                                'computing GA', 'computing EmA',
                                'computing KA'])

    out, err = check_output('python3 sfepy/examples/homogenization/rs_correctors.py -n')
    eok += report(out, '...', -2, -1, '1.644e-01', match_numbers=True)

    out, err = check_output('python3 sfepy/examples/large_deformation/compare_elastic_materials.py -n')
    eok += report(out, '...', -8, 5, '1.068759e-14', eps=1e-13)

    out, err = check_output('python3 sfepy/examples/linear_elasticity/linear_elastic_interactive.py')
    eok += report(out, '...', -18, 0, '1.62128841139e-14', eps=1e-13)

    out, err = check_output('python3 sfepy/examples/linear_elasticity/elastodynamic_identification.py --opt-conf=xtol=0.5')
    eok += report(out, '...', -5, 0, '2', match_numbers=True)

    out, err = check_output('python3 sfepy/examples/linear_elasticity/modal_analysis.py')
    eok += report(out, '...', -12, 5, '12142.11470773', eps=1e-13)

    out, err = check_output('python3 sfepy/examples/multi_physics/thermal_electric.py')
    eok += report(out, '...', -9, 5, '2.612933e-14', eps=1e-13)

    out, err = check_output('python3 sfepy/examples/diffusion/laplace_refine_interactive.py output')
    eok += report(out, '...', -3, 5, '2.675866e-15', eps=1e-13)

    out, err = check_output('python3 sfepy/examples/diffusion/laplace_iga_interactive.py -o output')
    eok += report(out, '...', -3, 5, '1.028134e-13', eps=1e-12)

    out, err = check_output('python3 sfepy/examples/dg/imperative_burgers_1D.py -o output')
    eok += report(out, '...', -3, 3, 'moment_1D_limiter')

    out, err = check_output('mpiexec -n 2 python3 sfepy/examples/diffusion/poisson_parallel_interactive.py output-parallel -2 --silent -ksp_monitor -options_left 0')
    eok += report(out, '...', -2, 4, '8.021313824020e-07', eps=1e-6)

    out, err = check_output('mpiexec -n 2 python3 sfepy/examples/multi_physics/biot_parallel_interactive.py output-parallel -2 --silent -ksp_monitor -options_left 0')
    eok += report(out, '...', -2, 4, '3.787214380277e-09', eps=1e-7)

    t1 = time.time()

    out, err = check_output("python3 -c \"import sfepy; sfepy.test('-v', '--disable-warnings')\"")
    tok, failed = report_tests(out, return_item=True)
    tok = {True : 'ok', False : 'fail'}[tok]

    t2 = time.time()

    fd = open('test_install_times.log', 'a+')
    fd.write('%s: examples: %.2f [s] (%d), tests: %.2f [s] (%s: %s)\n'
             % (time.ctime(t0), t1 - t0, eok, t2 - t1, tok, failed))
    fd.close()

if __name__ == '__main__':
    main()

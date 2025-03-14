import os.path as op
import inspect

import numpy as nm
import numpy.linalg as nla

from sfepy.base.base import Struct, assert_, IndexedStruct

def report(*argc):
    """All tests should print via this function."""
    format = '...' + ' %s' * len(argc)
    msg = format % argc
    print(msg)

def eval_coor_expression(expression, coor):

    x = coor[:, 0]
    y = coor[:, 1]
    if coor.shape[1] == 3:
        z = coor[:, 2]
    else:
        z = None

    env = {'x' : x, 'y' : y, 'z' : z}
    out = eval(expression, nm.__dict__, env)

    if isinstance(out, float):
        aux = nm.empty(coor.shape[0], dtype=nm.float64)
        aux.fill(out)
        out = aux

    return out

def compare_vectors(vec1, vec2, allowed_error=1e-8,
                    label1='vec1', label2='vec2', norm=None):

    diff_norm = nla.norm(vec1 - vec2, ord=norm)
    report('||%s - %s||: %e' % (label1, label2, diff_norm))
    if diff_norm > allowed_error:
        return False
    else:
        return True

def assert_equal(a, b, msg='assertion of equality failed!'):
    import scipy.sparse

    assert_base_types = (int, float, str, bytes, complex,
                         None.__class__, type, nm.number)
    if a is b: return

    def assert_dict(a, b):
        assert_(set(a.keys()) == set(b.keys()), msg)
        for i in a:
            assert_equal(a[i], b[i], msg)

    def assert_list(a, b):
        assert_(len(a) == len(b), msg)
        for i, j in zip(a, b):
            assert_equal(i, j)

    assert_(a.__class__ is b.__class__, msg)
    if isinstance(a, (int, float, str, assert_base_types, bytes, complex)):
        assert_(a == b, msg)

    elif isinstance(a, dict):
        assert_dict(a, b)

    elif isinstance(a, (list, tuple)):
        assert_list(a, b)

    elif isinstance(a, nm.ndarray):
        nm.testing.assert_array_equal(a,b)

    elif isinstance(a, (scipy.sparse.csr_matrix, scipy.sparse.csc_matrix)):
        nm.testing.assert_array_equal(a.data, b.data)
        nm.testing.assert_array_equal(a.indices, b.indices)
        nm.testing.assert_array_equal(a.indptr, b.indptr)

    elif isinstance(a, object):
        cls = a.__class__

        if hasattr(cls, '__slots__'):
            ad = dict((i,getattr(a,i)) for i in cls.__slots__)
            bd = dict((i,getattr(b,i)) for i in cls.__slots__)

        elif hasattr(a, '__dict__'):
            ad = a.__dict__
            bd = b.__dict__

        else:
            def members(obj):
                out = inspect.getmembers(obj, lambda x: not
                                         inspect.isroutine(x) )
                out = dict( (k,v) for k,v in out if not k.startswith('__'))
                return out
            ad = members(a)
            bd = members(b)

        assert_dict(ad, bd)

class NLSStatus(IndexedStruct):
    """
    Custom nonlinear solver status storing stopping condition of all
    time steps.
    """
    def __setitem__(self, key, val):
        IndexedStruct.__setitem__(self, key, val)
        if key == 'condition':
            self.conditions.append(val)

def check_conditions(conditions):
    ok = (conditions == 0).all()
    if not ok:
        report('nls stopping conditions:')
        report(conditions)
    return ok

def run_declaratice_example(ex_filename, define_args=None,
                            output_dir='', ext='.vtk',
                            remove_prefix=''):
    """
    Run a declarative example in `ex_filename` given relatively to
    ``sfepy.base_dir``.
    """
    import sfepy
    from sfepy.applications import solve_pde

    report('solving %s...' % ex_filename)

    if remove_prefix and ex_filename.startswith(remove_prefix):
        output_name = ex_filename.replace(remove_prefix, '')

    else:
        output_name = ex_filename

    output_name = op.splitext(output_name.replace('/', '-'))[0]
    filename = op.join(sfepy.base_dir, ex_filename)
    name = op.splitext(op.split(output_name)[1])[0]
    fmt = ext.replace('.', '')
    options = Struct(output_filename_trunk=name,
                     output_format=fmt if fmt != '' else 'vtk',
                     save_ebc=False, save_ebc_nodes=False,
                     save_regions=False,
                     save_regions_as_groups=False,
                     solve_not=False)
    status = IndexedStruct(nls_status=NLSStatus(conditions=[]))

    solve_pde(filename, define_args=define_args, options=options, status=status,
              output_dir=output_dir)
    report('%s solved' % ex_filename)

    return nm.array(status.nls_status.conditions)

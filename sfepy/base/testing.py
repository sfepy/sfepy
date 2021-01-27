from __future__ import print_function
import sys
import inspect

import numpy as nm
import numpy.linalg as nla

from sfepy.base.base import Struct, output, assert_

class TestCommon(Struct):

    @staticmethod
    def xfail(test_method):
        """
        Decorator that allows a test to fail.
        """
        def wrapper(self):
            try:
                ret = test_method(self)
            except:
                if self.debug:
                    raise
                ret = False

            if not ret:
                print('--- test expected to fail.')
                ret = True

            return ret

        return wrapper

    def get_number(self):
        methods = inspect.getmembers(self, inspect.ismethod)
        tests = [ii for ii in methods
                 if (len(ii[0]) > 5) and ii[0][:5] == 'test_']
        return len(tests)

    def run(self, debug=False, ifile=None):
        self.debug = debug

        ok = True
        n_fail = 0

        methods = inspect.getmembers(self, inspect.ismethod)
        if hasattr(self, 'tests'):
            dmethods = {}
            for key, method in methods:
                dmethods[key] = method
            tests = [(ii, dmethods[ii]) for ii in self.tests]
            print(tests)
        else:
            tests = [ii for ii in methods
                     if (len(ii[0]) > 5) and ii[0][:5] == 'test_']

        orig_prefix = output.get_output_prefix()
        for itest, (test_name, test_method) in enumerate(tests):
            if len(tests) > 1 and ifile is not None:
                output.set_output_prefix('[%d/%d] %s'\
                    % (ifile, itest + 1, orig_prefix))
            else:
                output.set_output_prefix('[%d] %s' % (ifile, orig_prefix))

            aux = ' %s: ' % test_name

            try:
                ret = test_method()
            except:
                if debug:
                    raise
                ret = False

            if not ret:
                aux = '---' + aux + 'failed!'
                n_fail += 1
                ok = False
            else:
                aux = '+++' + aux + 'ok'

            print(aux)

        output.set_output_prefix(orig_prefix)

        return ok, n_fail, len(tests)

    @staticmethod
    def report(*argc):
        """All tests should print via this function."""
        format = '...' + ' %s' * len(argc)
        msg = format % argc
        print(msg)

    @staticmethod
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

    @staticmethod
    def compare_vectors(vec1, vec2, allowed_error=1e-8,
                        label1='vec1', label2='vec2', norm=None):

        diff_norm = nla.norm(vec1 - vec2, ord=norm)
        TestCommon.report('||%s - %s||: %e' % (label1, label2, diff_norm))
        if diff_norm > allowed_error:
            return False
        else:
            return True

    @staticmethod
    def assert_equal(a, b, msg='assertion of equality failed!'):
        import scipy.sparse

        assert_base_types = (int, float, str, bytes, complex,
                             None.__class__, type)
        if sys.version_info < (3, 0):
            assert_base_types = assert_base_types + (unicode,) # NOQA
            #
            # Dirty fix for Py27/win64/NumPy combo (LK) NOQA
            #
            assert_base_types = assert_base_types + (long,)

        if a is b: return

        def assert_dict(a, b):
            assert_(set(a.keys()) == set(b.keys()), msg)
            for i in a:
                TestCommon.assert_equal(a[i], b[i], msg)

        def assert_list(a, b):
            assert_(len(a) == len(b), msg)
            for i, j in zip(a, b):
                TestCommon.assert_equal(i, j)

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

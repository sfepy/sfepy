from __future__ import absolute_import
from sfepy.base.testing import TestCommon
from sfepy.homogenization.engine import HomogenizationEngine as he
from sfepy.homogenization.engine import HomogenizationWorkerMulti as hwm
import six
import numpy as nm

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        test = Test(conf=conf, options=options)
        return test

    def test_dependencies(self):
        get_deps = hwm.get_sorted_dependencies

        coefs = {'A' : {'requires' : ['a', 'd', 'c.B']},
                 'B' : {'requires' : ['b']}}
        requirements = {'a' : {'requires' : ['b', 'c']},
                        'b' : {'requires' : ['c']},
                        'c' : {},
                        'd' : {'requires' : ['b', 'a']}}

        deps = get_deps(requirements, coefs, None)
        ok = ((deps == ['c', 'b', 'a', 'd', 'c.B', 'c.A'])
              or (deps == ['c', 'b', 'c.B', 'a', 'd', 'c.A'])
              or (deps == ['c', 'b', 'a', 'c.B', 'd', 'c.A']))
        self.report(deps, ':', ok)

        coefs['B']['requires'] = ['b', 'c.A']

        try:
            deps = get_deps(requirements, coefs, None)

        except ValueError as err:
            self.report('detected:', str(err))
            _ok = 'circular requirement "c.' in str(err)

        else:
            _ok = False

        self.report('circular dependency detection 1:', _ok)
        ok = ok and _ok

        coefs['B']['requires'] = ['b']
        requirements['c']['requires'] = ['d']

        try:
            deps = get_deps(requirements, coefs, None)

        except ValueError as err:
            self.report('detected:', str(err))
            _ok = 'circular requirement' in str(err)

        else:
            _ok = False

        self.report('circular dependency detection 2:', _ok)
        ok = ok and _ok

        return ok

    def test_chunk_micro(self):
        coefs = {'A' : {'requires' : ['a', 'd', 'c.B']},
                 'B' : {'requires' : ['b']}}
        requirements = {'a' : {'requires' : ['b', 'c']},
                        'b' : {'requires' : ['c']},
                        'c' : {},
                        'd' : {'requires' : ['b', 'a']}}

        volumes = {'total': {'expression': ''}}
        coefs = he.define_volume_coef(coefs, volumes)
        orig_deps_num = len(requirements) + len(coefs)

        num_workers, num_micro, chunks_per_worker = 5, 61, 2
        store_micro_idxs = [0, 1, 18, 20, 21]
        micro_chunk_tab, requirements, coefs = \
            hwm.chunk_micro_coors(num_workers, num_micro, requirements, coefs,
                                  chunks_per_worker, store_micro_idxs)

        dep_names = hwm.get_sorted_dependencies(requirements, coefs, None)

        ok = (orig_deps_num * len(micro_chunk_tab)) == len(dep_names)
        self.report('splitting into chunks:', ok)

        deps = {}
        for k in dep_names:
            chunk_id = int(k[-3:])
            nmic = len(range(*micro_chunk_tab[chunk_id].indices(num_micro)))
            deps[k] = [1] * nmic
            if k[2:] in coefs and 'Volume_total' not in k:
                reqs = '#'.join(coefs[k[2:]]['requires'])
                ok = ok and 'Volume_total' in reqs

        self.report('volume dependecy:', ok)

        deps = hwm.dechunk_reqs_coefs(deps, len(micro_chunk_tab))

        ok = ok and\
            nm.all([(nm.sum(v) == num_micro) for v in six.itervalues(deps)])
        self.report('merging chunks:', ok)

        return ok

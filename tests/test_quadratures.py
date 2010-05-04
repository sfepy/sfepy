import numpy as nm

from sfepy.base.testing import TestCommon
from sfepy.base.base import ordered_iteritems

class Test(TestCommon):

    @staticmethod
    def from_conf( conf, options ):

        test = Test(conf=conf, options=options)
        return test

    def test_weight_consistency(self):
        """
        Test if integral of 1 (= sum of weights) gives the domain volume.
        """
        from sfepy.fem.quadratures import quadrature_tables

        ok = True
        for geometry, qps in ordered_iteritems(quadrature_tables):
            self.report('geometry:', geometry)
            for order, qp in ordered_iteritems(qps):
                _ok = nm.allclose(qp.weights.sum(), qp.volume,
                                  rtol=0.0, atol=1e-15)
                self.report('order %d: %s' % (order, _ok))

                ok = ok and _ok

        return ok

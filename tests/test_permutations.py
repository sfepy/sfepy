import os.path as op
from sfepy.base.testing import TestCommon, assert_

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        return Test(conf = conf, options = options)

    def test_rcm(self):
        from sfepy.linalg import rcm, permute_in_place, save_sparse_txt
        from sfepy.fem import Mesh

        filename = '../database/tests/triquad.mesh'

        self.report('testing reversed Cuthill-McKee permutation')

        conf_dir = op.dirname(__file__)
        mesh = Mesh.from_file(filename, prefix_dir=conf_dir)

        graph = mesh.create_conn_graph()
        graph0 = graph.copy()

        save_sparse_txt(op.join(self.options.out_dir, 'test_rcm_graph_orig'),
                        graph, fmt='%d %d %d\n')

        perm = rcm(graph)

        permute_in_place(graph, perm)
        save_sparse_txt(op.join(self.options.out_dir, 'test_rcm_graph_rcm'),
                        graph, fmt='%d %d %d\n')
        
        assert_((graph0.indptr != graph.indptr).any())
        assert_((graph0.indices != graph.indices).any())

        permute_in_place(graph, perm, inverse=True)
        save_sparse_txt(op.join(self.options.out_dir, 'test_rcm_graph_rcm_inv'),
                        graph, fmt='%d %d %d\n')

        assert_((graph0.indptr == graph.indptr).all())
        assert_((graph0.indices == graph.indices).all())

        return True

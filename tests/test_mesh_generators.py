import os.path as op

from sfepy.base.testing import TestCommon

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        test = Test(conf=conf, options=options)
        return test

    def test_gen_block_mesh(self):
        from sfepy.mesh.mesh_generators import gen_block_mesh

        mesh = gen_block_mesh([1, 2, 3], [4, 3, 5], [2, 1, 0], verbose=False)
        filename = op.join(self.options.out_dir, 'gen_block.mesh')
        mesh.write(filename)

        self.report('block mesh generated')
        return True

    def test_gen_cylinder_mesh(self):
        from sfepy.mesh.mesh_generators import gen_cylinder_mesh

        mesh = gen_cylinder_mesh([0.5, 1, 2, 1.5, 3], [5, 4, 3], [0, 2, 1],
                                 axis='z', non_uniform=True, verbose=False)
        filename = op.join(self.options.out_dir, 'gen_cylinger.mesh')
        mesh.write(filename)

        self.report('cylinder mesh generated')
        return True

    def test_gen_extended_block_mesh(self):
        from sfepy.mesh.mesh_generators import gen_extended_block_mesh

        mesh = gen_extended_block_mesh([2, 3, 1], [5, 3, 4], [14, 10, 20],
                                       5, [1, 0, 2])
        filename = op.join(self.options.out_dir, 'gen_extended_block.mesh')
        mesh.write(filename)

        self.report('extended block mesh generated')
        return True

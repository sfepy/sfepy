from sfepy.mesh.mesh_generators import Mish_mesh_creator, gen_block_mesh
from sfepy.base.base import output

class Cube(Mish_mesh_creator):
      arg = Mish_mesh_creator.arg
      arg_formats =  (int, float, int, float, int, float)
      arg_defaults = (20, 5, arg(0), arg(1), arg(0), arg(1))
      filename_pattern = 'cube-%i_%.2f-%i_%.2f-%i_%.2f'

      def __call__(self, filename):
        output('creating new cube mesh')
        output('(%i nodes in %.2f) x (%i nodes in %.2f) x (%i nodes in %.2f)'
               % self.args)
        output('to file %s...' % filename)

        mesh = gen_block_mesh(self.args[1::2], self.args[0::2],
                              (0.0, 0.0, 0.0), name=filename)
        mesh.write(filename, io='auto')
        output('...done')

from sfepy.mesh.mesh_generators import Mish_mesh_creator
from sfepy.mesh.excube import excube

class Excube(Mish_mesh_creator):
      arg = Mish_mesh_creator.arg
      arg_formats =  (float, float, float, float, float, float, float)
      arg_defaults = (0.1, 3, 3, arg(1), arg(2), arg(1), arg(2))
      filename_pattern = 'excube-%.3f-%.2f+%.2f-%.2f+%.2f-%.2f+%.2f'

      def __call__(self, filename):
          excube(self.args[0], self.args[1::2], self.args[2::2], filename)
 				  

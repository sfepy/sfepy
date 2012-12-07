from sfepy.mesh.mesh_generators import Mish_mesh_creator
from sfepy.base.base import output
import subprocess
from sfepy import data_dir

class Sphere(Mish_mesh_creator):
      arg = Mish_mesh_creator.arg
      arg_formats =  (float, int, float)
      arg_defaults = (5, 41, arg(0))
      filename_pattern = 'cube-%i_%.2f-%i_%.2f-%i_%.2f' 


      def force_suffix(self):
        return '.mesh' 

      def __call__(self, filename):
        output('creating new sphere mesh (%i nodes, r=%.2f) and gradation %d'
               % args)
        output('to file %s...' % filename)
        defdir = os.path.join(data_dir, 'meshes')

        f = open(os.path.join(defdir, 'quantum', 'sphere.geo'))
        tmpfile = os.path.join(data_dir, 'tmp', 'sphere.geo.temp')
        ff = open(tmpfile, "w")
        ff.write("""
R = %i.0;
n = %i.0;
dens = %f;
""" % args)
        ff.write(f.read())
        f.close()
        ff.close()
        subprocess.call(['gmsh', '-3', tmpfile, '-format', 'mesh',
                         '-o', filename])

        output('...done')


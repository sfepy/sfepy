import os.path as op
import shutil

import numpy as nm

from sfepy.base.base import output, get_default, Struct
from sfepy.homogenization.coefficients import Coefficients
from sfepy.homogenization.coefs_base import MiniAppBase
from sfepy.homogenization.engine import HomogenizationEngine
from sfepy.applications import SimpleApp

class Volume(MiniAppBase):

    def __call__(self, problem=None):
        problem = get_default(problem, self.problem)
        problem.select_variables(self.variables)

        volume = problem.evaluate(self.expression)

        return volume

class HomogenizationApp( HomogenizationEngine ):

    def process_options( options ):
        """Application options setup. Sets default values for missing
        non-compulsory options."""
        get = options.get_default_attr
        
        print_digits = get( 'print_digits', 3 )

        float_format = get( 'float_format', '%8.3e' )
        coefs_filename = get( 'coefs_filename', 'coefs' )
        tex_names = get( 'tex_names', None )
        
        coefs = get( 'coefs', None, 'missing "coefs" in options!' )
        requirements = get( 'requirements', None,
                            'missing "requirements" in options!' )

        volume = get('volume', None)
        volumes = get('volumes', None)
        if volume is None and volumes is None:
            raise ValueError('missing "volume" in options!')

        return Struct( **locals() )
    process_options = staticmethod( process_options )

    def __init__(self, conf, options, output_prefix, **kwargs):
        SimpleApp.__init__(self, conf, options, output_prefix,
                           init_equations=False)

        self.setup_options()
        self.cached_coefs = None

        output_dir = self.problem.output_dir

        if conf._filename is not None:
            shutil.copyfile(conf._filename,
                            op.join(output_dir, op.basename(conf._filename)))

    def setup_options( self ):
        SimpleApp.setup_options( self )
        po = HomogenizationApp.process_options
        self.app_options += po( self.conf.options )

    def call( self, ret_all = False ):
        opts = self.app_options

        volume = {}
        if opts.volumes is not None:
            for vk, vv in opts.volumes.iteritems():
                if 'value' in vv:
                    volume[vk] = nm.float64(vv['value'])
                else:
                    volume[vk] = Volume('volume', self.problem, vv)()
        else:
            if opts.volume is not None:
                if 'value' in opts.volume:
                    vol = nm.float64(opts.volume['value'])
                else:
                    vol = Volume('volume', self.problem, opts.volume)()
                volume['total'] = vol

        for vk, vv in volume.iteritems():
            output('volume: %s = %.2f' % (vk, vv))

        if hasattr(opts.options, 'return_all'):
            ret_all = opts.options.return_all

        he = HomogenizationEngine( self.problem, self.options, volume = volume )

        aux = he( ret_all = ret_all)
        if ret_all:
            coefs, dependencies = aux
        else:
            coefs = aux

        coefs = Coefficients( **coefs.to_dict() )
        coefs.volume = volume
        
        prec = nm.get_printoptions()[ 'precision']
        if hasattr( opts, 'print_digits' ):
            nm.set_printoptions( precision = opts.print_digits )
        print coefs
        nm.set_printoptions( precision = prec )
##        pause()

        coef_save_name = op.join( opts.output_dir, opts.coefs_filename )
        coefs.to_file_hdf5( coef_save_name + '.h5' )
        coefs.to_file_txt( coef_save_name + '.txt',
                           opts.tex_names,
                           opts.float_format )

        if ret_all:
            return coefs, dependencies
        else:
            return coefs

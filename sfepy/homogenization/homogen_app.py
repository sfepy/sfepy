import os.path as op
import shutil

import numpy as nm

from sfepy.base.base import output, get_default, Struct
from sfepy.homogenization.coefficients import Coefficients
from sfepy.homogenization.coefs_base import MiniAppBase
from sfepy.homogenization.engine import HomogenizationEngine
from sfepy.applications import PDESolverApp

class Volume(MiniAppBase):

    def __call__(self, problem=None):
        problem = get_default(problem, self.problem)
        problem.select_variables(self.variables)

        volume = problem.evaluate(self.expression)

        return volume

def get_volume_from_options(options, problem):
    volume = {}

    if hasattr(options, 'volumes') and (options.volumes is not None):
        for vk, vv in options.volumes.iteritems():
            if 'value' in vv:
                volume[vk] = nm.float64(vv['value'])
            else:
                volume[vk] = Volume('volume', problem, vv)()

    elif hasattr(options, 'volume') and (options.volume is not None):
            if 'value' in options.volume:
                vol = nm.float64(options.volume['value'])
            else:
                vol = Volume('volume', problem, options.volume)()
            volume['total'] = vol

    return volume

class HomogenizationApp( HomogenizationEngine ):

    @staticmethod
    def process_options(options):
        """
        Application options setup. Sets default values for missing
        non-compulsory options.
        """
        get = options.get_default_attr

        volume = get('volume', None)
        volumes = get('volumes', None)
        if volume is None and volumes is None:
            raise ValueError('missing "volume" in options!')

        return Struct(print_digits=get('print_digits', 3),
                      float_format=get('float_format', '%8.3e'),
                      coefs_filename=get('coefs_filename', 'coefs'),
                      tex_names=get('tex_names', None),
                      coefs=get('coefs', None, 'missing "coefs" in options!'),
                      requirements=get('requirements', None,
                                       'missing "requirements" in options!'),
                      return_all=get('return_all', False),
                      volume=volume,
                      volumes=volumes)

    def __init__(self, conf, options, output_prefix, **kwargs):
        PDESolverApp.__init__(self, conf, options, output_prefix,
                              init_equations=False)

        self.setup_options()
        self.cached_coefs = None

        output_dir = self.problem.output_dir

        if conf._filename is not None:
            shutil.copyfile(conf._filename,
                            op.join(output_dir, op.basename(conf._filename)))

    def setup_options( self ):
        PDESolverApp.setup_options(self)
        po = HomogenizationApp.process_options
        self.app_options += po( self.conf.options )

    def call(self, ret_all=False, verbose=False):
        opts = self.app_options

        volume = get_volume_from_options(opts, self.problem)

        for vk, vv in volume.iteritems():
            output('volume: %s = %.2f' % (vk, vv))

        he = HomogenizationEngine( self.problem, self.options, volume = volume )

        aux = he( ret_all = ret_all)
        if ret_all:
            coefs, dependencies = aux
        else:
            coefs = aux

        coefs = Coefficients( **coefs.to_dict() )
        coefs.volume = volume

        if verbose:
            prec = nm.get_printoptions()[ 'precision']
            if hasattr(opts, 'print_digits'):
                nm.set_printoptions(precision=opts.print_digits)
            print coefs
            nm.set_printoptions(precision=prec)

        coef_save_name = op.join( opts.output_dir, opts.coefs_filename )
        coefs.to_file_hdf5( coef_save_name + '.h5' )
        coefs.to_file_txt( coef_save_name + '.txt',
                           opts.tex_names,
                           opts.float_format )

        if ret_all:
            return coefs, dependencies
        else:
            return coefs

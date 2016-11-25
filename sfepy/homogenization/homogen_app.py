from __future__ import print_function
from __future__ import absolute_import
import os.path as op
import shutil

import numpy as nm

from sfepy.base.base import get_default, Struct
from sfepy.homogenization.coefficients import Coefficients
from sfepy.homogenization.engine import HomogenizationEngine
from sfepy.applications import PDESolverApp
import sfepy.base.multiproc as multiproc
import sfepy.discrete.fem.periodic as per
import sfepy.linalg as la
from six.moves import range

class Volume(MiniAppBase):

    def __call__(self, problem=None):
        problem = get_default(problem, self.problem)
        problem.select_variables(self.variables)

        volume = problem.evaluate(self.expression)

        return volume

def get_volume_from_options(options, problem, ncoors=None):
    def get_vol(opt_vol, problem):
        if ncoors is None:
            if 'value' in opt_vol:
                out = nm.float64(opt_vol['value'])
            else:
                out = Volume('volume', problem, opt_vol)()
        else:
            n = ncoors.shape[0]
            out = nm.empty((n,), dtype=nm.float64)
            for im in range(n):
                if 'value' in opt_vol:
                    out[im] = nm.float64(opt_vol['value'])
                else:
                    problem.set_mesh_coors(ncoors[im], update_fields=True,
                                           clear_all=False, actual=True)
                    out[im] = Volume('volume', problem, opt_vol)()

        return out

    volume = {}
    if hasattr(options, 'volumes') and (options.volumes is not None):
        for vk, vv in six.iteritems(options.volumes):
            volume[vk] = get_vol(vv, problem)

    elif hasattr(options, 'volume') and (options.volume is not None):
        volume['total'] = get_vol(options.volume, problem)

    return volume

class HomogenizationApp(HomogenizationEngine):

    @staticmethod
    def process_options(options):
        """
        Application options setup. Sets default values for missing
        non-compulsory options.
        """
        get = options.get

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
                      macro_deformation=get('macro_deformation', None),
                      mesh_update_variable=get('mesh_update_variable', None),
                      mesh_update_corrector=get('mesh_update_corrector', None),
                      time_tag=get('time_tag', ''),
                      volume=volume,
                      volumes=volumes)

    def __init__(self, conf, options, output_prefix, **kwargs):
        PDESolverApp.__init__(self, conf, options, output_prefix,
                              init_equations=False)

        self.setup_options()
        self.cached_coefs = None
        self.n_micro = kwargs.get('n_micro', None)
        self.macro_deformation = None
        self.micro_coors = None
        self.updating_corrs = None

        mac_def = self.app_options.macro_deformation
        if mac_def is not None and isinstance(mac_def, nm.ndarray):
            self.n_micro = mac_def.shape[0]
            self.setup_macro_deformation(mac_def)

        if self.n_micro is not None:
            coors = self.problem.domain.get_mesh_coors()
            self.micro_coors = nm.empty((self.n_micro,) + coors.shape,
                                        dtype=nm.float64)
            for im in range(self.n_micro):
                self.micro_coors[im,...] = coors

        output_dir = self.problem.output_dir

        if conf._filename is not None:
            shutil.copyfile(conf._filename,
                            op.join(output_dir, op.basename(conf._filename)))

    def setup_options(self):
        PDESolverApp.setup_options(self)
        po = HomogenizationApp.process_options
        self.app_options += po(self.conf.options)

    def setup_macro_deformation(self, mtx_F):
        """
        Setup macroscopic deformation gradient.
        """
        self.macro_deformation = mtx_F

    def update_micro_coors(self, ret_val=False):
        """
        Update microstructures coordinates according to the deformation
        gradient and corrector functions.
        """
        dim = self.macro_deformation.shape[1]
        mtx_e = self.macro_deformation - nm.eye(dim)
        ncoors = self.micro_coors
        ncoors += la.dot_sequences(ncoors, mtx_e, 'ABT')
        if self.updating_corrs is not None:
            upd_var = self.app_options.mesh_update_variable
            for ii, corr in enumerate(self.updating_corrs):
                update_corr = nm.array(\
                    [corr.states[jj][upd_var] for jj in corr.components]).T
                gg = mtx_e[ii,...].reshape((dim**2, 1))
                ncoors[ii] += nm.dot(update_corr, gg).reshape(ncoors[ii].shape)

        if ret_val:
            return ncoors

    def call(self, verbose=False, ret_all=None, time_tag=''):
        """
        Call the homogenization engine and compute the homogenized
        coefficients.

        Parameters
        ----------
        verbose : bool
            If True, print the computed coefficients.
        ret_all : bool or None
            If not None, it can be used to override the 'return_all' option.
            If True, also the dependencies are returned.
        time_tag: str
            The time tag used in file names.

        Returns
        -------
        coefs : Coefficients instance
            The homogenized coefficients.
        dependencies : dict
            The dependencies, if `ret_all` is True.
        """
        opts = self.app_options

        ret_all = get_default(ret_all, opts.return_all)

        if not hasattr(self, 'he'):
            self.he = HomogenizationEngine(self.problem, self.options)

        if self.micro_coors is not None:
            self.he.set_micro_coors(self.update_micro_coors(ret_val=True))

        volume = get_volume_from_options(opts, self.problem, self.micro_coors)
        self.he.set_volume(volume)
        for vk, vv in six.iteritems(volume):
            output('volume: %s = ' % vk, vv)

        if multiproc.use_multiprocessing:
            upd_var = self.app_options.mesh_update_variable
            if upd_var is not None:
                uvar = self.problem.create_variables([upd_var])[upd_var]
                uvar.field.mappings0 = multiproc.get_dict('mappings0')
            per.periodic_cache = multiproc.get_dict('periodic_cache')

        aux = self.he(ret_all=ret_all, time_tag=time_tag)
        if ret_all:
            coefs, dependencies = aux
            # store correctors for coors update
            if opts.mesh_update_corrector is not None:
                self.updating_corrs =\
                    dependencies[opts.mesh_update_corrector]
        else:
            coefs = aux

        coefs = Coefficients(**coefs.to_dict())
        coefs.volume = volume

        if verbose:
            prec = nm.get_printoptions()[ 'precision']
            if hasattr(opts, 'print_digits'):
                nm.set_printoptions(precision=opts.print_digits)
            print(coefs)
            nm.set_printoptions(precision=prec)

        coef_save_name = op.join(opts.output_dir, opts.coefs_filename)
        coefs.to_file_hdf5(coef_save_name + '%s.h5' % time_tag)
        coefs.to_file_txt(coef_save_name + '%s.txt' % time_tag,
                          opts.tex_names,
                          opts.float_format)

        if ret_all:
            return coefs, dependencies
        else:
            return coefs

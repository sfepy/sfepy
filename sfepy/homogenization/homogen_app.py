import os.path as op
import shutil

import numpy as nm

from sfepy.base.base import get_default, Struct
from sfepy.homogenization.coefficients import Coefficients
from sfepy.homogenization.engine import HomogenizationEngine
from sfepy.applications import PDESolverApp
import sfepy.discrete.fem.periodic as per
import sfepy.linalg as la
import sfepy.base.multiproc as multi


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
                      mesh_update_variable=get('mesh_update_variable', None),
                      macro_data=get('macro_data', None),
                      micro_update=get('micro_update', {}),
                      n_micro=get('n_micro', None),
                      multiprocessing=get('multiprocessing', True),
                      use_mpi=get('use_mpi', False),
                      store_micro_idxs=get('store_micro_idxs', []),
                      volume=volume,
                      volumes=volumes)

    def __init__(self, conf, options, output_prefix, **kwargs):
        PDESolverApp.__init__(self, conf, options, output_prefix,
                              init_equations=False)

        self.setup_options()
        self.n_micro = kwargs.get('n_micro',
                                  self.app_options.get('n_micro', None))
        self.updating_corrs = None
        self.micro_state_cache = {}
        self.multiproc_mode = None
        self.micro_states = None if self.n_micro is None else {}

        # macroscopic data given in problem options dict.
        macro_data = self.app_options.macro_data
        if macro_data is not None:
            self.n_micro = macro_data[list(macro_data.keys())[0]].shape[0]
            self.setup_macro_data(macro_data)

        if self.n_micro is not None:
            for k in self.app_options.micro_update:
               if not k == 'coors':
                   self.micro_states[k] = None

            coors = self.problem.domain.get_mesh_coors()
            c_sh = (self.n_micro,) + coors.shape
            self.micro_states['coors'] = nm.empty(c_sh, dtype=nm.float64)

            mac_ids = kwargs.get('mac_ids',
                                  self.app_options.get('mac_ids', None))
            self.micro_states['id'] = []
            for im in range(self.n_micro):
                self.micro_states['coors'][im] = coors
                self.micro_states['id'].append(
                    mac_ids[im] if mac_ids is not None else im)

        output_dir = self.problem.output_dir

        if conf._filename is not None:
            shutil.copyfile(conf._filename,
                            op.join(output_dir, op.basename(conf._filename)))

    def setup_options(self):
        PDESolverApp.setup_options(self)
        po = HomogenizationApp.process_options
        self.app_options += po(self.conf.options)
        if hasattr(self, 'he'):
            self.he.setup_options()

    def setup_macro_data(self, data):
        """
        Setup macroscopic deformation gradient.
        """
        self.macro_data = data
        self.problem.homogenization_macro_data = self.macro_data

    def get_micro_cache_key(self, key, icoor, itime):
        tt = '' if itime is None else '_t%03d' % itime
        return '%s_%d%s' % (key, icoor, tt)

    def update_micro_states(self):
        """
        Update microstructures state according to the macroscopic data
        and corrector functions.
        """
        def calculate_local_update(state, corrs, var, macro_vals, mul=1):
            for ic, corr in enumerate(corrs):
                if state is None:
                    sh = corr.states[corr.components[0]][var].shape \
                        if hasattr(corr, 'states') else corr.state[var].shape
                    state = nm.zeros((len(corrs),) + sh, dtype=nm.float64)
                else:
                    sh = state[ic].shape

                if hasattr(corr, 'states'):
                    corr_arr = nm.array(
                        [corr.states[jj][var] for jj in corr.components]).T
                    mval = macro_vals[ic].reshape((corr_arr.shape[1], 1))
                    state[ic] += mul * nm.dot(corr_arr, mval).reshape(sh)
                else:
                    if macro_vals is None:
                        state[ic] += mul * corr.state[var].reshape(sh)
                    else:
                        state[ic] += mul\
                            * (corr.state[var] * macro_vals[ic]).reshape(sh)

            return state

        micro_update = self.app_options.micro_update
        for key, upd_obj in micro_update.items():
            if '_prev' in key:
                continue

            state = self.micro_states[key]

            if key + '_prev' in micro_update and state is not None:
                self.micro_states[key + '_prev'] = state.copy()

            if key == 'coors':
                if hasattr(upd_obj, '__call__'):
                    upd_obj(state, self.macro_data, self.problem)
                else:
                    # macro strain - in the first sequence of the list
                    mtx_e = self.macro_data[upd_obj[0][2]]
                    state += la.dot_sequences(state, mtx_e, 'ABT')

            if hasattr(upd_obj, '__call__'):
                upd_obj(state, self.macro_data, self.problem)
            else:
                if self.updating_corrs is not None:
                    for v in upd_obj:
                        if len(v) == 4:
                            cname, vname, mname, mul = v
                        else:
                            cname, vname, mname = v
                            mul = 1

                        macro_data = None if mname is None \
                            else self.macro_data[mname]
                        if cname is not None:
                            state0 = calculate_local_update(
                                state, self.updating_corrs[cname],
                                vname, macro_data, mul)
                        else:
                            state += macro_data[..., 0]

                        if state0 is not state:
                            self.micro_states[key] = state0
                            state = state0

    def call(self, verbose=False, ret_all=None, itime=None, iiter=None):
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

        force_init_he = hasattr(self.problem, 'force_init_he')\
            and self.problem.force_init_he
        if not hasattr(self, 'he') or force_init_he:
            volumes = {}
            if hasattr(opts, 'volumes') and (opts.volumes is not None):
                volumes.update(opts.volumes)
            elif hasattr(opts, 'volume') and (opts.volume is not None):
                volumes['total'] = opts.volume
            else:
                volumes['total'] = 1.0

            self.he = HomogenizationEngine(self.problem, self.options,
                                           volumes=volumes)

        if self.micro_states is not None:
            self.update_micro_states()
            self.he.set_micro_states(self.micro_states)

        multiproc_mode = None
        if opts.multiprocessing and multi.use_multiprocessing:
            multiproc, multiproc_mode = multi.get_multiproc(mpi=opts.use_mpi)

            if multiproc_mode is not None:
                upd_var = self.app_options.mesh_update_variable
                if upd_var is not None:
                    uvar = self.problem.create_variables([upd_var])[upd_var]
                    uvar.field.mappings0 = multiproc.get_dict('mappings0',
                                                              soft_set=True)
                per.periodic_cache = multiproc.get_dict('periodic_cache',
                                                        soft_set=True)

        time_tag = ('' if itime is None else '_t%03d' % itime)\
            + ('' if iiter is None else '_i%03d' % iiter)

        aux = self.he(ret_all=ret_all, time_tag=time_tag)
        if ret_all:
            coefs, dependencies = aux
            # store correctors for coors update
            self.updating_corrs = {}
            for upd_obj in opts.micro_update.values():
                if upd_obj is not None and not hasattr(upd_obj, '__call__'):
                    for v in upd_obj:
                        cr = v[0]
                        if cr is not None:
                            self.updating_corrs[cr] = dependencies[cr]
        else:
            coefs = aux

        if coefs is not None:
            coefs = Coefficients(**coefs.to_dict())

            if verbose:
                prec = nm.get_printoptions()['precision']
                if hasattr(opts, 'print_digits'):
                    nm.set_printoptions(precision=opts.print_digits)
                print(coefs)
                nm.set_printoptions(precision=prec)

            ms_cache = self.micro_state_cache
            for ii in self.app_options.store_micro_idxs:
                for k in self.micro_states.keys():
                    key = self.get_micro_cache_key(k, ii, itime)
                    ms_cache[key] = self.micro_states[k][ii]

            coef_save_name = op.join(opts.output_dir, opts.coefs_filename)
            coefs.to_file_hdf5(coef_save_name + '%s.h5' % time_tag)
            coefs.to_file_txt(coef_save_name + '%s.txt' % time_tag,
                              opts.tex_names,
                              opts.float_format)

        if ret_all:
            return coefs, dependencies
        else:
            return coefs

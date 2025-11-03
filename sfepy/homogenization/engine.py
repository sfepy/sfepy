from copy import copy
from functools import partial
import os.path as osp
from sfepy.base.base import output, get_default, Struct
from sfepy.applications import PDESolverApp, Application
from .coefs_base import MiniAppBase, CoefEval
from sfepy.discrete.evaluate import eval_equations
import sfepy.homogenization.multiproc as multiproc


def insert_sub_reqs(reqs, levels, req_info):
    """Recursively build all requirements in correct order."""
    all_reqs = []
    for _, req in enumerate(reqs):
        # Coefficients are referenced as 'c.<name>'...
        areq = req[2:] if req.startswith('c.') else req

        try:
            rargs = req_info[areq]
        except KeyError:
            raise ValueError('requirement "%s" is not defined!' % req)

        sub_reqs = rargs.get('requires', [])

        if req in levels:
            raise ValueError('circular requirement "%s"!' % (req))

        if sub_reqs:
            levels.append(req)
            sreqs = insert_sub_reqs(sub_reqs, levels, req_info)
            all_reqs += [ii for ii in sreqs if ii not in all_reqs]
            levels.pop()

        if req not in all_reqs:
            all_reqs.append(req)

    return all_reqs


def get_dict_idxval(dict_array, idx):
    return {k: v[idx] for k, v in dict_array.items()}


def deps_to_corrs(problem, dependencies, keys=None):
    if keys is None:
        keys = dependencies.keys()

    vs = problem.get_variables(auto_create=True)
    out = {}
    for k in keys:
        if not k.startswith('c.'):
            v = dependencies[k]
            if hasattr(v, 'get_output'):
                out[k] = v.get_output(variables=vs)
            else:
                out[k] = v

    return out


def rm_chtag(key):
    idx = key.find('|ch:')
    return key[:idx] if idx > 0 else key


def chunk_req_info(name, req_info):
    idx = name.find('|ch:')
    if idx > 0:
        name0, ctag = name[:idx], name[idx:]
        out = req_info[name0].copy()
        reqs = out.get('requires', [])
        if len(reqs) > 0:
            out['requires'] = [k + ctag for k in reqs]

        return out
    else:
        return req_info[name]


class CoefVolume(MiniAppBase):
    def __call__(self, volume, problem=None, data=None):
        problem = get_default(problem, self.problem)

        term_mode = self.term_mode
        equations, variables = problem.create_evaluable(self.expression,
                                                        term_mode=term_mode)
        return eval_equations(equations, variables, term_mode=term_mode)


class HomogenizationWorker:
    def __call__(self, problem, options, post_process_hook,
                 req_info, coef_info, micro_states, store_micro_idxs,
                 chunk_size=None, time_tag=''):
        """Calculate homogenized correctors and coefficients.

        Parameters
        ----------
        problem : problem
            The problem definition - microscopic problem.
        opts : struct
            The options of the homogenization application.
        post_process_hook : function
            The postprocessing hook.
        req_info : dict
            The definition of correctors.
        coef_info : dict
            The definition of homogenized coefficients.
        micro_states : array
            The configurations of multiple microstructures.
        store_micro_idxs : list of int
            The indices of microstructures whose results are to be stored.
        chunk_size: int or None
            The number of microstructures in one chunks
        time_tag : str
            The label corresponding to the actual time step and iteration,
            used in the corrector file names.

        Returns
        -------
        dependencies : dict
            The computed correctors and coefficients.
        save_names : list
            The names of computed dependencies.
        """
        dependencies = {}
        save_names = {}

        sorted_names = self.get_sorted_dependencies(req_info, coef_info,
                                                    options.compute_only)
        for name in sorted_names:
            if not name.startswith('c.'):
                if micro_states is not None:
                    req_info[name]['store_idxs'] = (store_micro_idxs, 0)

            val, sn = self.calculate_req(problem, options, post_process_hook,
                                         req_info, coef_info, dependencies,
                                         micro_states, time_tag,
                                         None, name, '0')

            dependencies[name] = val
            save_names.update(sn)

        return dependencies, save_names

    @staticmethod
    def get_sorted_dependencies(req_info, coef_info, compute_only):
        "Make corrs and coefs list sorted according to the dependencies."

        reqcoef_info = copy(coef_info)
        reqcoef_info.update(req_info)
        compute_names = set(get_default(compute_only, list(coef_info.keys())))
        compute_names = ['c.' + key for key in compute_names]

        dep_names = []
        for coef_name in compute_names:
            requires = coef_info[coef_name[2:]].get('requires', [])
            deps = insert_sub_reqs(copy(requires), [], reqcoef_info)\
                + [coef_name]
            for dep in deps:
                if dep not in dep_names:
                    dep_names.append(dep)

        return dep_names

    @staticmethod
    def calculate(mini_app, problem, dependencies, dep_requires,
                  save_names, micro_states, chunk_tab, mode, proc_id):
        data = {rm_chtag(key): dependencies[key] for key in dep_requires
                if 'Volume_' not in key}
        volume = {rm_chtag(key[9:]): dependencies[key] for key in dep_requires
                  if 'Volume_' in key}
        mini_app.requires = [rm_chtag(ii) for ii in mini_app.requires
                             if 'c.Volume_' not in ii]

        if micro_states is None:
            if mode == 'coefs':
                val = mini_app(volume, data=data)
            else:
                if mini_app.save_name is not None:
                    save_names[mini_app.name] = mini_app.get_save_name_base()
                val = mini_app(data=data)
        else:
            if '|ch:' in mini_app.name and chunk_tab is not None:
                chunk_id = int(mini_app.name.split('|ch:')[1])
                chunk_tag = f'-{chunk_id}'
                local_state = \
                    {k: v[chunk_tab[chunk_id]] if v is not None else None
                     for k, v in micro_states.items()}
            else:
                chunk_tag = ''
                local_state = micro_states

            val = []
            if hasattr(mini_app, 'store_idxs') and mode == 'reqs':
                save_name = mini_app.save_name

            local_coors = local_state['coors']
            for im in range(len(local_coors)):
                output(f'== micro {proc_id}{chunk_tag}-{im} ==')
                problem.micro_state = (local_state, im)
                problem.set_mesh_coors(local_coors[im], update_fields=True,
                                       clear_all=False, actual=True)

                if mode == 'coefs':
                    val.append(mini_app(get_dict_idxval(volume, im),
                                        data=get_dict_idxval(data, im)))
                else:
                    if hasattr(mini_app, 'store_idxs')\
                            and im in mini_app.store_idxs[0]:
                        store_id = '_%04d' % (mini_app.store_idxs[1] + im)
                        if save_name is not None:
                            mini_app.save_name = save_name + store_id
                            key = mini_app.name
                            if key in save_names:
                                save_names[key].append(
                                    mini_app.get_save_name_base())
                            else:
                                save_names[key] =\
                                    [mini_app.get_save_name_base()]
                    else:
                        mini_app.save_name = None

                    val.append(mini_app(data=get_dict_idxval(data, im)))

                    if len(val) == 1 and val[0].name == 'update_coors':
                        local_coors[im] += val[0].state

        return val

    @staticmethod
    def calculate_req(problem, opts, post_process_hook,
                      req_info, coef_info, dependencies,
                      micro_states, time_tag, chunk_tab, name, proc_id):
        """Calculate a requirement, i.e. correctors or coefficients.

        Parameters
        ----------
        problem : problem
            The problem definition related to the microstructure.
        opts : struct
            The options of the homogenization application.
        post_process_hook : function
            The postprocessing hook.
        name : str
            The name of the requirement.
        req_info : dict
            The definition of correctors.
        coef_info : dict
            The definition of homogenized coefficients.
        dependencies : dict
            The dependencies required by the correctors/coefficients.
        micro_states : array
            The configurations of multiple microstructures.
        time_tag : str
            The label corresponding to the actual time step and iteration,
            used in the corrector file names.
        chunk_tab : list
            In the case of multiprocessing the requirements are divided into
            several chunks that are solved in parallel.
        proc_id : int
            The id number of the processor (core) which is solving the actual
            chunk.

        Returns
        -------
        val : coefficient/corrector or list of coefficients/correctors
            The resulting homogenized coefficients or correctors.
        save_names : dict
            The dictionary containing names of saved correctors.
        """
        save_names = {}

        # compute coefficient
        if name.startswith('c.'):
            output('computing %s...' % name[2:])

            cargs = chunk_req_info(name[2:], coef_info)
            mini_app = MiniAppBase.any_from_conf(name[2:], problem, cargs)
            problem.clear_equations()

            # Pass only the direct dependencies, not the indirect ones.
            dep_requires = cargs.get('requires', [])
            val = HomogenizationWorker.calculate(mini_app, problem,
                                                 dependencies, dep_requires,
                                                 save_names, micro_states,
                                                 chunk_tab, 'coefs', proc_id)

            output('...done')

        # compute corrector(s)
        else:
            output('computing dependency %s...' % name)

            rargs = chunk_req_info(name, req_info)
            mini_app = MiniAppBase.any_from_conf(name, problem, rargs)
            mini_app.setup_output(save_formats=opts.save_formats,
                                  post_process_hook=post_process_hook,
                                  split_results_by=opts.split_results_by)
            if mini_app.save_name is not None:
                mini_app.save_name += time_tag

            problem.clear_equations()

            # Pass only the direct dependencies, not the indirect ones.
            dep_requires = rargs.get('requires', [])
            val = HomogenizationWorker.calculate(mini_app, problem,
                                                 dependencies, dep_requires,
                                                 save_names, micro_states,
                                                 chunk_tab, 'reqs', proc_id)

            output('...done')

        return val, save_names

    @staticmethod
    def recover_micro(problem, corrs, rhook, macro_data):
        out = []
        corrs_dict = deps_to_corrs(problem, corrs)
        for local_macro, label in macro_data:
            output(label, verbose=True)
            out.append(rhook(problem, corrs_dict, local_macro))

        return out

    @staticmethod
    def process_reqs_coefs(orig, nchunk, store_idxs=[]):
        new = {}
        for k, v in orig.items():
            if k == 'filenames':
                continue

            for ii in range(nchunk):
                lab = f'|ch:{ii}'
                key = k + lab
                val = new[key] = {}
                if 'requires' in v:
                    val['requires'] = [jj + lab for jj in v['requires']]
                if len(store_idxs) > 0:
                    if len(store_idxs[ii][0]) > 0:
                        val['store_idxs'] = store_idxs[ii]
                    else:
                        val['save_name'] = None

        return new

class HomogenizationWorkerMulti(HomogenizationWorker):
    @staticmethod
    def calculate_req(*args):
        proc_id = ''.join(k for k in multiproc.get_proc_id() if k.isdigit())
        output.set_output_prefix(f'he-w{proc_id}:')
        new_args = args + (proc_id, )

        return HomogenizationWorker.calculate_req(*new_args)

    @staticmethod
    def recover_req(*args):
        proc_id = ''.join(k for k in multiproc.get_proc_id() if k.isdigit())
        output.set_output_prefix(f'he-w{proc_id}:')
        new_args = args + (proc_id, )

        return HomogenizationWorker.calculate_req(*new_args)

    @staticmethod
    def rhook_call(problem, dependencies, rhook, macro_data):
        proc_id = ''.join(k for k in multiproc.get_proc_id() if k.isdigit())
        output.set_output_prefix(f'he-w{proc_id}:')

        corrs = deps_to_corrs(problem, dependencies)
        local_macro, label = macro_data
        output(label, verbose=True)

        return rhook(problem, corrs, local_macro)

    @staticmethod
    def recover_micro(problem, corrs, rhook, macro_data):
        dependencies = multiproc.get_dict('dependencies')
        dependencies.update({k: v for k, v in corrs.items()
                             if k not in dependencies})

        workers = multiproc.get_workers()
        if workers is None:
            conf_file_dir = osp.split(problem.conf.__file__)[0]
            max_workers = problem.conf.options.get('max_workers', None)
            workers = multiproc.init_workers(conf_file_dir, max_workers)

        rhook = problem.conf.options.get('recovery_hook', None)
        rhook = problem.conf.get_function(rhook)
        wfun = partial(HomogenizationWorkerMulti.rhook_call,
                       problem, dependencies, rhook)

        out = [k for k in workers.map(wfun, macro_data)]

        return out

    def __call__(self, problem, options, post_process_hook,
                 req_info, coef_info, micro_states, store_micro_idxs,
                 chunk_size, time_tag=''):
        """Calculate homogenized correctors and coefficients in separated
           processes.

        Parameters
        ----------
        The same parameters as :class:`HomogenizationWorker`, extended by:
        chunk_size : int
            The number of chunks per one worker.

        Returns
        -------
        The same returns as :class:`HomogenizationWorker`.
        """
        dependencies = multiproc.get_dict('dependencies', clear=True)

        if micro_states is not None and chunk_size is not None:
            micro_chunk_tab, req_info_ch, coef_info_ch = \
                self.chunk_micro_tasks(len(micro_states['coors']),
                                       req_info, coef_info,
                                       chunk_size, store_micro_idxs)
        else:
            micro_chunk_tab = None
            req_info_ch, coef_info_ch = req_info, coef_info

        workers = multiproc.get_workers()
        if workers is None:
            conf_file_dir = osp.split(problem.conf.__file__)[0]
            max_workers = problem.conf.options.get('max_workers', None)
            workers = multiproc.init_workers(conf_file_dir, max_workers)

        wfun = partial(HomogenizationWorkerMulti.calculate_req, problem,
                       options, post_process_hook, req_info, coef_info,
                       dependencies, micro_states, time_tag, micro_chunk_tab)

        sorted_names = self.get_sorted_dependencies(req_info_ch, coef_info_ch,
                                                    options.compute_only)

        # calculate number of dependencies and inverse map
        numdeps = {}
        inverse_deps = {}
        for name in sorted_names:
            if name.startswith('c.'):
                reqs = coef_info_ch[name[2:]].get('requires', [])
            else:
                reqs = req_info_ch[name].get('requires', [])
            numdeps[name] = len(reqs)
            if len(reqs) > 0:
                for req in reqs:
                    if req in inverse_deps:
                        inverse_deps[req].append(name)
                    else:
                        inverse_deps[req] = [name]

        save_names = {}
        remaining = len(sorted_names)

        while remaining > 0:
            if len(numdeps) > 0:
                tasks = [k for k, v in numdeps.items() if v == 0]
                numdeps = {k: v for k, v in numdeps.items() if v > 0}

            for task, (dep, snames) in zip(tasks, workers.map(wfun, tasks)):
                dependencies[task] = dep
                save_names.update(snames)

                if task in inverse_deps:
                    for itask in inverse_deps[task]:
                        numdeps[itask] -= 1  # itask depends on task

                remaining -= 1

        if micro_chunk_tab is not None:
            dependencies = self.dechunk_reqs_coefs(dependencies,
                                                   len(micro_chunk_tab))
        else:
            dependencies = dependencies.copy()

        return dependencies, save_names

    @staticmethod
    def chunk_micro_tasks(num_micro, reqs, coefs,
                          chunk_size, store_micro_idxs=[]):
        """
        Split multiple microproblems into several chunks
        that can be processed in parallel.

        Parameters
        ----------
        num_micro : int
            The number of microstructures.
        reqs : dict
            The requirement definitions.
        coefs : dict
            The coefficient definitions.
        chunk_size : int
            The number of microproblems per chunk.
        store_micro_idxs : list of int
            The indices of microstructures whose results are to be stored.

        Returns
        -------
        micro_tab : list of slices
            The indices of microproblems contained in each chunk.
        new_reqs : dict
            The new requirement definitions - .
        new_coefs : dict
            The new coefficient definitions.
        """
        micro_tab = []
        store_idxs = []
        for ii in range(0, num_micro, chunk_size):
            jj = chunk_size + ii
            micro_tab.append(slice(ii, min([num_micro, jj])))
            if len(store_micro_idxs) > 0:
                store_idxs.append(([k - ii for k in store_micro_idxs
                                    if k >= ii and k < jj], ii))

        nchunk = len(micro_tab)
        self = HomogenizationWorkerMulti
        new_reqs = self.process_reqs_coefs(reqs, nchunk, store_idxs)
        new_coefs = self.process_reqs_coefs(coefs, nchunk)

        return micro_tab, new_reqs, new_coefs

    @staticmethod
    def dechunk_reqs_coefs(deps, num_chunks):
        """
        Merge the results related to the multiple microproblems.

        Parameters
        ----------
        deps : dict
            The calculated dependencies.
        num_chunks : int
            The number of chunks.

        Returns
        -------
        new_deps : dict
            The merged dependencies.
        """
        new_deps = {}
        for dep in set([rm_chtag(k) for k in deps.keys()]):
            new_deps[dep] = sum([deps[f'{dep}|ch:{ii}']
                                 for ii in range(num_chunks)], [])

        return new_deps


class HomogenizationEngine(PDESolverApp):
    @staticmethod
    def process_options(options):
        get = options.get

        return Struct(coefs=get('coefs', None,
                                'missing "coefs" in options!'),
                      requirements=get('requirements', None,
                                       'missing "requirements" in options!'),
                      compute_only=get('compute_only', None),
                      multiprocessing=get('multiprocessing', True),
                      store_micro_idxs=get('store_micro_idxs', []),
                      chunk_size=get('chunk_size', 1),
                      chunks_per_worker=get('chunks_per_worker', None),
                      save_formats=get('save_formats', ['vtk', 'h5']),
                      coefs_info=get('coefs_info', None))

    def __init__(self, problem, options, app_options=None,
                 volumes=None, output_prefix='he:', **kwargs):
        """Bypasses PDESolverApp.__init__()!"""
        Application.__init__(self, problem.conf, options, output_prefix,
                             **kwargs)
        self.problem = problem
        self.setup_options(app_options=app_options)
        self.setup_output_info(self.problem, self.options)
        self.volumes = volumes
        self.micro_states = None

    def setup_options(self, app_options=None):
        PDESolverApp.setup_options(self)
        app_options = get_default(app_options, self.conf.options)

        po = HomogenizationEngine.process_options
        self.app_options += po(app_options)

    def set_micro_states(self, states):
        self.micro_states = states

    @staticmethod
    def define_volume_coef(coef_info, volumes):
        """
        Define volume coefficients and make all other dependent on them.

        Parameters
        ----------
        coef_info : dict
            The coefficient definitions.
        volumes : dict
            The definitions of volumes.

        Returns
        -------
        coef_info : dict
            The coefficient definitions extended by the volume coefficients.
        """
        vcfkeys = []
        cf_vols = {}
        for vk, vv in volumes.items():
            cfkey = 'Volume_%s' % vk
            vcfkeys.append('c.' + cfkey)
            if 'value' in vv:
                cf_vols[cfkey] = {'expression': '%e' % float(vv['value']),
                                  'class': CoefEval}
            else:
                cf_vols[cfkey] = {'expression': vv['expression'],
                                  'class': CoefVolume}

        for cf in coef_info.values():
            if 'requires' in cf:
                cf['requires'] += vcfkeys
            else:
                cf['requires'] = vcfkeys

        coef_info.update(cf_vols)

        return coef_info

    def recover(self, corrs, rhook, macro_data):
        if self.app_options.multiprocessing and multiproc.use_multiprocessing:
            rcall = HomogenizationWorkerMulti.recover_micro
        else:
            rcall = HomogenizationWorker.recover_micro

        return rcall(self.problem, corrs, rhook, macro_data)

    def call(self, ret_all=False, time_tag=''):
        problem = self.problem
        opts = self.app_options

        # Some coefficients can require other coefficients - resolve their
        # order here.
        req_info = getattr(self.conf, opts.requirements, {})
        coef_info = getattr(self.conf, opts.coefs, {})
        coef_info = self.define_volume_coef(coef_info, self.volumes)

        is_store_filenames = coef_info.pop('filenames', None) is not None
        if opts.multiprocessing and multiproc.use_multiprocessing:
            worker = HomogenizationWorkerMulti()
            if (self.micro_states is not None
                    and opts.chunks_per_worker is not None):
                nch = (multiproc.cpu_count() - 1) * opts.chunks_per_worker
                # ceil division
                chunk_size = -(len(self.micro_states['coors']) // -nch)
            else:
                chunk_size = opts.chunk_size
        else:  # no multiprocessing
            worker = HomogenizationWorker()
            chunk_size = None

        dependencies, save_names = worker(problem, opts,
                                          self.post_process_hook,
                                          req_info, coef_info,
                                          self.micro_states,
                                          opts.store_micro_idxs,
                                          chunk_size, time_tag)

        deps = {}

        if dependencies is not None:
            coefs = Struct()
            for name in dependencies.keys():
                data = dependencies[name]
                if name.startswith('c.'):
                    coef_name = name[2:]
                    cstat = coef_info[coef_name].get('status', 'main')
                    # remove "auxiliary" coefs
                    if not cstat == 'auxiliary':
                        setattr(coefs, coef_name, data)
                else:
                    deps[name] = data

            # Store filenames of all requirements as a "coefficient".
            if is_store_filenames and save_names is not None:
                coefs.save_names = save_names

            if opts.coefs_info is not None:
                coefs.info = opts.coefs_info

        if ret_all:
            return coefs, deps
        else:
            return coefs

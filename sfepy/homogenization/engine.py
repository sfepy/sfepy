import gc
from copy import copy

from sfepy.base.base import output, get_default, Struct
from sfepy.applications import PDESolverApp, Application
from .coefs_base import MiniAppBase, CoefEval
from .utils import rm_multi
from sfepy.discrete.evaluate import eval_equations
import sfepy.base.multiproc as multi
import numpy as nm


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


class CoefVolume(MiniAppBase):
    def __call__(self, volume, problem=None, data=None):
        problem = get_default(problem, self.problem)

        term_mode = self.term_mode
        equations, variables = problem.create_evaluable(self.expression,
                                                        term_mode=term_mode)
        return eval_equations(equations, variables, term_mode=term_mode)


class HomogenizationWorker:
    def __call__(self, problem, options, post_process_hook,
                 req_info, coef_info,
                 micro_states, store_micro_idxs, time_tag=''):
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

            val = self.calculate_req(problem, options, post_process_hook,
                                     name, req_info, coef_info, save_names,
                                     dependencies, micro_states,
                                     time_tag)

            dependencies[name] = val
            gc.collect()

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
        if micro_states is None:
            data = {key: dependencies[key] for key in dep_requires
                    if 'Volume_' not in key}
            volume = {key[9:]: dependencies[key]
                      for key in dep_requires if 'Volume_' in key}
            mini_app.requires = [ii for ii in mini_app.requires
                                 if 'c.Volume_' not in ii]

            if mode == 'coefs':
                val = mini_app(volume, data=data)
            else:
                if mini_app.save_name is not None:
                    save_names[mini_app.name] = mini_app.get_save_name_base()
                val = mini_app(data=data)
        else:
            data = {rm_multi(key): dependencies[key]
                    for key in dep_requires if 'Volume_' not in key}
            volume = {rm_multi(key[9:]): dependencies[key]
                      for key in dep_requires if 'Volume_' in key}
            mini_app.requires = [ii for ii in mini_app.requires
                                 if 'c.Volume_' not in ii]

            if '|multiprocessing_' in mini_app.name\
                    and chunk_tab is not None:
                chunk_id = int(mini_app.name[-3:])
                chunk_tag = '-%d' % (chunk_id + 1)
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
                output('== micro %s%s-%d =='
                       % (proc_id, chunk_tag, im + 1))
                problem.micro_state = (local_state, im)
                problem.set_mesh_coors(local_coors[im], update_fields=True,
                                       clear_all=False, actual=True)

                if mode == 'coefs':
                    val.append(mini_app(get_dict_idxval(volume, im),
                                        data=get_dict_idxval(data, im)))
                else:
                    if hasattr(mini_app, 'store_idxs')\
                            and im in mini_app.store_idxs[0]:
                        store_id = '_%04d'\
                            % (mini_app.store_idxs[1] + im)
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
                      name, req_info, coef_info, save_names, dependencies,
                      micro_states, time_tag='', chunk_tab=None, proc_id='0'):
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
        save_names : dict
            The dictionary containing names of saved correctors.
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
        """
        # compute coefficient
        if name.startswith('c.'):
            coef_name = name[2:]

            output('computing %s...' % coef_name)

            cargs = coef_info[coef_name]
            mini_app = MiniAppBase.any_from_conf(coef_name, problem, cargs)
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

            rargs = req_info[name]
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

        return val


class HomogenizationWorkerMulti(HomogenizationWorker):
    def __init__(self, num_workers):
        self.num_workers = num_workers

    def __call__(self, problem, options, post_process_hook,
                 req_info, coef_info,
                 micro_states, store_micro_idxs, chunks_per_worker,
                 time_tag=''):
        """Calculate homogenized correctors and coefficients.

        Parameters
        ----------
        The same parameters as :class:`HomogenizationWorker`, extended by:
        chunks_per_worker : int
            The number of chunks per one worker.

        Returns
        -------
        The same returns as :class:`HomogenizationWorker`.
        """
        multiproc = multi.multiproc_proc

        dependencies = multiproc.get_dict('dependecies', clear=True)
        save_names = multiproc.get_dict('save_names', clear=True)
        numdeps = multiproc.get_dict('numdeps', clear=True)
        remaining = multiproc.get_int_value('remaining', 0)
        tasks = multiproc.get_queue('tasks')
        lock = multiproc.get_lock('lock')

        if micro_states is not None:
            micro_chunk_tab, req_info, coef_info = \
                self.chunk_micro_tasks(self.num_workers,
                                       len(micro_states['coors']),
                                       req_info, coef_info,
                                       chunks_per_worker, store_micro_idxs)
        else:
            micro_chunk_tab = None

        sorted_names = self.get_sorted_dependencies(req_info, coef_info,
                                                    options.compute_only)

        remaining.value = len(sorted_names)

        # calculate number of dependencies and inverse map
        inverse_deps = {}
        for name in sorted_names:
            if name.startswith('c.'):
                reqs = coef_info[name[2:]].get('requires', [])
            else:
                reqs = req_info[name].get('requires', [])
            numdeps[name] = len(reqs)
            if len(reqs) > 0:
                for req in reqs:
                    if req in inverse_deps:
                        inverse_deps[req].append(name)
                    else:
                        inverse_deps[req] = [name]

        for name in sorted_names:
            if numdeps[name] == 0:
                tasks.put(name)

        workers = []
        for ii in range(self.num_workers):
            args = (tasks, lock, remaining, numdeps, inverse_deps,
                    problem, options, post_process_hook, req_info,
                    coef_info, save_names, dependencies, micro_states,
                    time_tag, micro_chunk_tab, str(ii + 1))
            w = multiproc.Process(target=self.calculate_req_multi,
                                  args=args)
            w.start()
            workers.append(w)

        # block until all workes are terminated
        for w in workers:
            w.join()

        if micro_states is not None:
            dependencies = self.dechunk_reqs_coefs(dependencies,
                                                   len(micro_chunk_tab))

        return dependencies, save_names

    @staticmethod
    def calculate_req_multi(tasks, lock, remaining, numdeps, inverse_deps,
                            problem, opts, post_process_hook,
                            req_info, coef_info, save_names, dependencies,
                            micro_states, time_tag, chunk_tab, proc_id):
        """Calculate a requirement in parallel.

        Parameters
        ----------
        tasks : queue
            The queue of requirements to be solved.
        lock : lock
            The multiprocessing lock used to ensure save access to the global
            variables.
        remaining : int
            The number of remaining requirements.
        numdeps : dict
            The number of dependencies for the each requirement.
        inverse_deps : dict
            The inverse dependencies - which requirements depend
            on a given one.

        For the definition of other parameters see 'calculate_req'.
        """
        while remaining.value > 0:
            name = tasks.get()

            if name is None:
                continue

            save_names_loc = {}
            val = HomogenizationWorker.calculate_req(problem, opts,
                post_process_hook, name, req_info, coef_info, save_names_loc,
                dependencies, micro_states, time_tag, chunk_tab, proc_id)

            lock.acquire()
            dependencies[name] = val
            remaining.value -= 1
            if name in inverse_deps:
                for iname in inverse_deps[name]:
                    numdeps[iname] -= 1  # iname depends on name
                    if numdeps[iname] == 0:  # computed all direct dependecies?
                        tasks.put(iname)  # yes, put iname to queue

            save_names.update(save_names_loc)
            lock.release()

    @staticmethod
    def process_reqs_coefs(old, num_workers, store_idxs=[]):
        new = {}
        for k, v in old.items():
            if k == 'filenames':
                new[k] = v.copy()
                continue

            for ii in range(num_workers):
                lab = '|multiprocessing_%03d' % ii
                key = k + lab
                new[key] = v.copy()
                val = new[key]
                if 'requires' in val:
                    val['requires'] = [jj + lab for jj in val['requires']]
                if len(store_idxs) > 0:
                    if len(store_idxs[ii][0]) > 0:
                        val['store_idxs'] = store_idxs[ii]
                    else:
                        val['save_name'] = None

        return new

    @staticmethod
    def chunk_micro_tasks(num_workers, num_micro, reqs, coefs,
                          chunks_per_worker=1, store_micro_idxs=[]):
        """
        Split multiple microproblems into several chunks
        that can be processed in parallel.

        Parameters
        ----------
        num_workers : int
            The number of available CPUs.
        num_micro : int
            The number of microstructures.
        reqs : dict
            The requirement definitions.
        coefs : dict
            The coefficient definitions.
        chunks_per_worker : int
            The number of chunks per one worker.
        store_micro_idxs : list of int
            The indices of microstructures whose results are to be stored.

        Returns
        -------
        micro_tab : list of slices
            The indices of microproblems contained in each chunk.
        new_reqs : dict
            The new requirement definitions.
        new_coefs : dict
            The new coefficient definitions.
        """
        chsize = int(nm.ceil(float(num_micro)
                     / (num_workers * chunks_per_worker)))

        micro_tab = []
        store_idxs = []
        for ii in range(0, num_micro, chsize):
            jj = chsize + ii
            chunk_end = num_micro if jj > num_micro else jj
            micro_tab.append(slice(ii, chunk_end))
            if len(store_micro_idxs) > 0:
                store_idxs.append(([k - ii for k in store_micro_idxs
                                    if k >= ii and k < jj], ii))

        nw = len(micro_tab)
        self = HomogenizationWorkerMulti
        new_reqs = self.process_reqs_coefs(reqs, nw, store_idxs)
        new_coefs = self.process_reqs_coefs(coefs, nw)

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
        for ii in range(num_chunks):
            ilab = '_%03d' % ii
            for k in deps.keys():
                idx = k.rfind('|multiprocessing_')
                if idx > 0:
                    if not(k[-4:] == ilab):
                        continue
                    key = k[:idx]
                    if key in new_deps:
                        new_deps[key] += deps[k]
                    else:
                        new_deps[key] = deps[k]
                else:
                    new_deps[k] = deps[k]

        return new_deps


class HomogenizationWorkerMultiMPI(HomogenizationWorkerMulti):
    def __call__(self, problem, options, post_process_hook,
                 req_info, coef_info,
                 micro_states, store_micro_idxs, chunks_per_worker,
                 time_tag=''):
        """Calculate homogenized correctors and coefficients.

        Parameters and Returns
        ----------------------
        The same parameters and returns as :class:`HomogenizationWorkerMulti`.
        """
        multiproc = multi.multiproc_mpi

        dependencies = multiproc.get_dict('dependecies', clear=True)
        save_names = multiproc.get_dict('save_names', clear=True)
        numdeps = multiproc.get_dict('numdeps', mutable=True, clear=True)
        remaining = multiproc.get_int_value('remaining', 0)
        tasks = multiproc.get_queue('tasks')

        if micro_states is not None:
            micro_chunk_tab, req_info, coef_info = \
                self.chunk_micro_tasks(self.num_workers,
                                       len(micro_states['coors']),
                                       req_info, coef_info,
                                       chunks_per_worker, store_micro_idxs)
        else:
            micro_chunk_tab = None

        sorted_names = self.get_sorted_dependencies(req_info, coef_info,
                                                    options.compute_only)

        # calculate number of dependencies and inverse map
        inverse_deps = {}
        loc_numdeps = {}
        for name in sorted_names:
            if name.startswith('c.'):
                reqs = coef_info[name[2:]].get('requires', [])
            else:
                reqs = req_info[name].get('requires', [])
            loc_numdeps[name] = len(reqs)
            if len(reqs) > 0:
                for req in reqs:
                    if req in inverse_deps:
                        inverse_deps[req].append(name)
                    else:
                        inverse_deps[req] = [name]

        if multiproc.mpi_rank == multiproc.mpi_master:  # master node
            for k, v in loc_numdeps.items():
                numdeps[k] = v

            remaining.value = len(sorted_names)

            for name in sorted_names:
                if numdeps[name] == 0:
                    tasks.put(name)

            multiproc.master_loop()
            multiproc.master_send_continue()

            if micro_states is not None:
                dependencies = self.dechunk_reqs_coefs(dependencies,
                                                       len(micro_chunk_tab))

            multiproc.master_send_task('deps', dependencies)
            multiproc.master_send_continue()

            return dependencies, save_names

        else:  # slave node
            lock = multiproc.RemoteLock()
            multiproc.slave_get_task('engine')

            self.calculate_req_multi(tasks, lock, remaining, numdeps,
                                     inverse_deps, problem, options,
                                     post_process_hook, req_info,
                                     coef_info, save_names, dependencies,
                                     micro_states,
                                     time_tag, micro_chunk_tab,
                                     str(multiproc.mpi_rank + 1))

            multiproc.slave_task_done('engine')
            multiproc.wait_for_tag(multiproc.tags.CONTINUE)
            task, deps = multiproc.slave_get_task('get_deps')
            multiproc.wait_for_tag(multiproc.tags.CONTINUE)

            return deps, None


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
                      use_mpi=get('use_mpi', False),
                      store_micro_idxs=get('store_micro_idxs', []),
                      chunks_per_worker=get('chunks_per_worker', 1),
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

    def call(self, ret_all=False, time_tag=''):
        problem = self.problem
        opts = self.app_options

        # Some coefficients can require other coefficients - resolve their
        # order here.
        req_info = getattr(self.conf, opts.requirements, {})
        coef_info = getattr(self.conf, opts.coefs, {})
        coef_info = self.define_volume_coef(coef_info, self.volumes)

        is_store_filenames = coef_info.pop('filenames', None) is not None

        multiproc_mode = None
        if opts.multiprocessing and multi.use_multiprocessing:
            multiproc, multiproc_mode = multi.get_multiproc(mpi=opts.use_mpi)
            if multiproc_mode == 'mpi':
                HomogWorkerMulti = HomogenizationWorkerMultiMPI
            elif multiproc_mode == 'proc':
                HomogWorkerMulti = HomogenizationWorkerMulti
            else:
                multiproc_mode = None

        if multiproc_mode is not None:
            num_workers = multi.get_num_workers()
            # if self.micro_states is not None:
            #     n_micro = len(self.micro_states['coors'])
            #     if num_workers > n_micro:
            #         num_workers = n_micro
            worker = HomogWorkerMulti(num_workers)
            dependencies, save_names = \
                worker(problem, opts, self.post_process_hook,
                       req_info, coef_info, self.micro_states,
                       self.app_options.store_micro_idxs,
                       self.app_options.chunks_per_worker, time_tag)

        else:  # no multiprocessing
            worker = HomogenizationWorker()
            dependencies, save_names = \
                worker(problem, opts, self.post_process_hook,
                       req_info, coef_info, self.micro_states,
                       self.app_options.store_micro_idxs, time_tag)

        deps = {}

        if save_names is None and dependencies is not None:  # slave mode
            coefs = None
            for name in dependencies.keys():
                data = dependencies[name]
                if not name.startswith('c.'):
                    deps[name] = data
        else:
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
            if is_store_filenames:
                for name in save_names.keys():
                    if '|multiprocessing_' in name:
                        mname = rm_multi(name)
                        if mname in save_names:
                            save_names[mname] += save_names[name]
                        else:
                            save_names[mname] = save_names[name]
                        del(save_names[name])

                if multiproc_mode == 'proc':
                    coefs.save_names = save_names._getvalue()
                else:
                    coefs.save_names = save_names


            if opts.coefs_info is not None:
                coefs.info = opts.coefs_info

        if ret_all:
            return coefs, deps
        else:
            return coefs

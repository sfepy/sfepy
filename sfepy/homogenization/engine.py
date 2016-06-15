from __future__ import absolute_import
from copy import copy

from sfepy.base.base import output, get_default, Struct
from sfepy.applications import PDESolverApp, Application
from .coefs_base import MiniAppBase

try:
    import multiprocessing
    import Queue
    _use_multiprocessing = multiprocessing.cpu_count() > 1
except:
    _use_multiprocessing = False

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

        if req in all_reqs:
            raise ValueError('circular requirement "%s"!' % (req))
        else:
            all_reqs.append(req)

    return all_reqs

class HomogenizationEngine(PDESolverApp):

    @staticmethod
    def process_options(options):
        get = options.get

        return Struct(coefs=get('coefs', None,
                                'missing "coefs" in options!'),
                      requirements=get('requirements', None,
                                       'missing "requirements" in options!'),
                      compute_only=get('compute_only', None),
                      save_format=get('save_format', 'vtk'),
                      dump_format=get('dump_format', 'h5'),
                      coefs_info=get('coefs_info', None))

    def __init__(self, problem, options, app_options=None,
                 volume=None, output_prefix='he:', **kwargs):
        """Bypasses PDESolverApp.__init__()!"""
        Application.__init__(self, problem.conf, options, output_prefix,
                             **kwargs)
        self.problem = problem
        self.setup_options(app_options=app_options)
        self.setup_output_info(self.problem, self.options)

        if volume is None:
            self.volume = self.problem.evaluate(self.app_options.total_volume)

        else:
            self.volume = volume

    def setup_options(self, app_options=None):
        PDESolverApp.setup_options(self)
        app_options = get_default(app_options, self.conf.options)

        po = HomogenizationEngine.process_options
        self.app_options += po(app_options)

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
    def calculate_req(problem, opts, volume, post_process_hook,
                      name, req_info, coef_info, sd_names, dependencies):
        # compute coefficient
        if name.startswith('c.'):
            coef_name = name[2:]

            output('computing %s...' % coef_name)

            cargs = coef_info[coef_name]
            mini_app = MiniAppBase.any_from_conf(coef_name, problem, cargs)
            problem.clear_equations()

            # Pass only the direct dependencies, not the indirect ones.
            dep_requires = cargs.get('requires', [])
            data = {key: dependencies[key] for key in dep_requires}

            val = mini_app(volume, data=data)

            output('...done')

        # compute corrector(s)
        else:
            output('computing dependency %s...' % name)

            rargs = req_info[name]
            mini_app = MiniAppBase.any_from_conf(name, problem, rargs)
            mini_app.setup_output(save_format=opts.save_format,
                                  dump_format=opts.dump_format,
                                  post_process_hook=post_process_hook,
                                  file_per_var=opts.file_per_var)

            if not '(not_set)' in mini_app.get_save_name_base():
                sd_names['s.' + mini_app.name] = mini_app.get_save_name_base()
            if not '(not_set)' in mini_app.get_dump_name_base():
                sd_names['d.' + mini_app.name] = mini_app.get_dump_name_base()

            problem.clear_equations()

            # Pass only the direct dependencies, not the indirect ones.
            dep_requires = rargs.get('requires', [])
            data = {key: dependencies[key] for key in dep_requires}

            val = mini_app(data=data)

            output('...done')

        return val

    @staticmethod
    def calculate_req_multi(tasks, lock, remaining, numdeps, inverse_deps,
                            problem, opts, volume, post_process_hook,
                            req_info, coef_info, sd_names, dependencies):
        while remaining.value > 0:
            try:
                name = tasks.get(False, 0.01) # get (or wait for) a task
            except Queue.Empty:
                break

            store = {}
            val = HomogenizationEngine.calculate_req(problem, opts,
                                                     volume, post_process_hook,
                                                     name, req_info, coef_info,
                                                     store, dependencies)
            lock.acquire()
            dependencies[name] = val
            remaining.value -= 1
            if name in inverse_deps:
                for ii, iname in inverse_deps[name]:
                    numdeps[ii] -= 1 # iname depends on name
                    if numdeps[ii] == 0: # computed all direct dependecies?
                        tasks.put(iname)  # yes, put iname to queue
            sd_names.update(store)
            lock.release()

    def call(self, ret_all=False):
        problem = self.problem
        opts = self.app_options

        # Some coefficients can require other coefficients - resolve their
        # order here.
        req_info = getattr(self.conf, opts.requirements, {})
        coef_info = getattr(self.conf, opts.coefs, {})

        is_store_filenames = coef_info.pop('filenames', None) is not None
        sorted_names = self.get_sorted_dependencies(req_info, coef_info,
                                                    opts.compute_only)

        use_multiprocessing = _use_multiprocessing\
            and getattr(self.conf.options, 'multiprocessing', True)\
            and len(sorted_names) > 2

        coefs = Struct()
        if use_multiprocessing:
            manager = multiprocessing.Manager()
            dependencies = manager.dict()
            sd_names = manager.dict()
            numdeps = manager.list()
            remaining = manager.Value('i', len(sorted_names))
            tasks = multiprocessing.Queue()
            lock = multiprocessing.Lock()

            # calculate namber of dependencies and inverse map
            inverse_deps = {}
            for ii, name in enumerate(sorted_names):
                if name.startswith('c.'):
                    reqs = coef_info[name[2:]].get('requires', [])
                else:
                    reqs = req_info[name].get('requires', [])
                numdeps.append(len(reqs))
                if len(reqs) > 0:
                    for req in reqs:
                        if req in inverse_deps:
                            inverse_deps[req].append((ii, name))
                        else:
                            inverse_deps[req] = [(ii, name)]

            for ii, name in enumerate(sorted_names):
                if numdeps[ii] == 0:
                    tasks.put(name)

            num_workers = multiprocessing.cpu_count()
            workers = []
            for ii in xrange(num_workers):
                args = (tasks, lock, remaining, numdeps, inverse_deps,
                        problem, opts, self.volume, self.post_process_hook,
                        req_info, coef_info, sd_names, dependencies)
                w = multiprocessing.Process(target=self.calculate_req_multi,
                                            args=args)
                w.start()
                workers.append(w)

            # block until all workes are terminated
            for w in workers:
                w.join()

        else: # no mlutiprocessing
            dependencies = {}
            sd_names = {}
            for name in sorted_names:
                val = self.calculate_req(problem, opts,
                                         self.volume, self.post_process_hook,
                                         name, req_info, coef_info,
                                         sd_names, dependencies)
                dependencies[name] = val

        coefs = Struct()
        deps = {}
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
            save_names = {}
            dump_names = {}
            for name in sd_names.keys():
                val = sd_names[name]
                if name.startswith('s.'):
                    save_names[name[2:]] = val
                elif name.startswith('d.'):
                    dump_names[name[2:]] = val
            coefs.save_names = save_names
            coefs.dump_names = dump_names

        if opts.coefs_info is not None:
            coefs.info = opts.coefs_info

        if ret_all:
            return coefs, deps
        else:
            return coefs

from copy import copy

from sfepy.base.base import output, get_default, Struct
from sfepy.applications import PDESolverApp, Application
from coefs_base import MiniAppBase

def insert_sub_reqs( reqs, levels, req_info ):
    """Recursively build all requirements in correct order."""
    all_reqs = []
##     print '>', levels, reqs
    for ii, req in enumerate( reqs ):
        # Coefficients are referenced as 'c.<name>'...
        areq = req
        if req.startswith('c.'):
            areq = req[2:]

        try:
            rargs = req_info[areq]
        except KeyError:
            raise ValueError('requirement "%s" is not defined!' % req)

        sub_reqs = rargs.get( 'requires', [] )
##         print '*', ii, req, sub_reqs

        if req in levels:
            raise ValueError('circular requirement "%s"!' % (req))

        if sub_reqs:
            levels.append( req )
            all_reqs.extend( insert_sub_reqs( sub_reqs, levels, req_info ) )
            levels.pop()
            
        if req in all_reqs:
            raise ValueError('circular requirement "%s"!' % (req))
        else:
            all_reqs.append( req )

##     print all_reqs
##     pause()
    return all_reqs

class HomogenizationEngine(PDESolverApp):

    @staticmethod
    def process_options(options):
        get = options.get_default_attr

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
        Application.__init__( self, problem.conf, options, output_prefix,
                              **kwargs )
        self.problem = problem
        self.setup_options(app_options=app_options)
        self.setup_output_info( self.problem, self.options )

        if volume is None:
            self.volume = self.problem.evaluate(self.app_options.total_volume)

        else:
            self.volume = volume

    def setup_options(self, app_options=None):
        PDESolverApp.setup_options(self)
        app_options = get_default(app_options, self.conf.options)

        po = HomogenizationEngine.process_options
        self.app_options += po(app_options)

    def compute_requirements( self, requirements, dependencies, store ):
        problem = self.problem

        opts = self.app_options
        req_info = getattr( self.conf, opts.requirements )

        requires = insert_sub_reqs( copy( requirements ), [], req_info )
        
        for req in requires:
            if req in dependencies and (dependencies[req] is not None):
                continue

            output( 'computing dependency %s...' % req )

            rargs = req_info[req]

            mini_app = MiniAppBase.any_from_conf( req, problem, rargs )
            mini_app.setup_output( save_format = opts.save_format,
                                   dump_format = opts.dump_format,
                                   post_process_hook = self.post_process_hook,
                                   file_per_var = opts.file_per_var )
            store( mini_app )

            problem.clear_equations()

            # Pass only the direct dependencies, not the indirect ones.
            dep_requires = rargs.get('requires', [])
            data = {}
            for key in dep_requires:
                data[key] = dependencies[key]

            dep = mini_app(data=data)

            dependencies[req] = dep
            output( '...done' )

        return dependencies
        
    def call( self, ret_all = False ):
        problem = self.problem

        opts = self.app_options
        coef_info = getattr( self.conf, opts.coefs )

        compute_names = set(get_default(opts.compute_only, coef_info.keys()))
        compute_names = ['c.' + key for key in compute_names]

        is_store_filenames = coef_info.pop('filenames', None) is not None
        try:
            compute_names.remove('c.filenames')
        except:
            pass

        dependencies = {}
        save_names = {}
        dump_names = {}
        def store_filenames( app ):
            if not '(not_set)' in app.get_save_name_base():
                save_names[app.name] = app.get_save_name_base()
            if not '(not_set)' in app.get_dump_name_base():
                dump_names[app.name] = app.get_dump_name_base()

        # Some coefficients can require other coefficients - resolve their
        # order here.
        req_info = self.conf.get(opts.requirements, {})
        info = copy(coef_info)
        info.update(req_info)
        all_deps = set(compute_names)
        sorted_names = []
        for coef_name in compute_names:
            cargs = coef_info[coef_name[2:]]
            requires = cargs.get('requires', [])
            deps = insert_sub_reqs(copy(requires), [], info)
            all_deps.update(deps)

            aux = [key for key in deps if key.startswith('c.')] + [coef_name]
            sorted_names.extend(aux)

        sorted_coef_names = []
        for name in sorted_names:
            if name[2:] not in sorted_coef_names:
                sorted_coef_names.append(name[2:])

        coefs = Struct()
        for coef_name in sorted_coef_names:
            cargs = coef_info[coef_name]
            output('computing %s...' % coef_name)
            requires = cargs.get('requires', [])
            requirements = [name for name in requires if not
                            name.startswith('c.')]

            self.compute_requirements(requirements, dependencies,
                                      store_filenames)

            for name in requires:
                if name.startswith('c.'):
                    dependencies[name] = getattr(coefs, name[2:])

            mini_app = MiniAppBase.any_from_conf( coef_name, problem, cargs )

            problem.clear_equations()

            # Pass only the direct dependencies, not the indirect ones.
            data = {}
            for key in requires:
                data[key] = dependencies[key]

            val = mini_app(self.volume, data=data)
            setattr(coefs, coef_name, val)
            output( '...done' )

        # remove "auxiliary" coefs
        for coef_name in sorted_coef_names:
            cstat = coef_info[coef_name].get('status', 'main')
            if cstat == 'auxiliary':
                delattr(coefs, coef_name)

        # Store filenames of all requirements as a "coefficient".
        if is_store_filenames:
            coefs.save_names = save_names
            coefs.dump_names = dump_names

        if opts.coefs_info is not None:
            coefs.info = opts.coefs_info

        if ret_all:
            return coefs, dependencies
        else:
            return coefs

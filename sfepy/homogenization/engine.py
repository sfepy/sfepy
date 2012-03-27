from copy import copy, deepcopy

from sfepy.base.base import output, get_default, Struct
from sfepy.applications import SimpleApp, Application
from sfepy.fem.region import sort_by_dependency
from coefs_base import MiniAppBase

def insert_sub_reqs( reqs, levels, req_info ):
    """Recursively build all requirements in correct order."""
    all_reqs = []
##     print '>', levels, reqs
    for ii, req in enumerate( reqs ):
        try:
            rargs = req_info[req]
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

class HomogenizationEngine( SimpleApp ):

    @staticmethod
    def process_options(options):
        get = options.get_default_attr

        return Struct(coefs=get('coefs', None,
                                'missing "coefs" in options!'),
                      requirements=get('requirements', None,
                                       'missing "requirements" in options!'),
                      save_format=get('save_format', 'vtk'),
                      dump_format=get('dump_format', 'h5'))

    def __init__(self, problem, options, app_options=None,
                 volume=None, output_prefix='he:', **kwargs):
        """Bypasses SimpleApp.__init__()!"""
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
        SimpleApp.setup_options(self)
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

        is_store_filenames = coef_info.pop('filenames', None) is not None

        dependencies = {}
        save_names = {}
        dump_names = {}
        def store_filenames( app ):
            if not '(not_set)' in app.get_save_name_base():
                save_names[app.name] = app.get_save_name_base()
            if not '(not_set)' in app.get_dump_name_base():
                dump_names[app.name] = app.get_dump_name_base()

        def _get_parents(req_list):
            out = []
            for req_name in req_list:
                aux = req_name.split('.')
                if len(aux) == 2:
                    out.append(aux[1])
            return out

        # Some coefficients can require other coefficients - resolve theirorder
        # here.
        graph = {}
        for coef_name, cargs in coef_info.iteritems():
            if not coef_name in graph:
                graph[coef_name] = [0]

            requires = cargs.get('requires', [])
            for parent in _get_parents(requires):
                graph[coef_name].append(parent)
                requires.remove('c.' + parent)
            
        sorted_coef_names = sort_by_dependency(deepcopy(graph))
        ## print graph
        ## print sorted_coef_names
        
        coefs = Struct()
        for coef_name in sorted_coef_names:
            cargs = coef_info[coef_name]
            output( 'computing %s...' % coef_name )
            requires = cargs.get( 'requires', [] )

            self.compute_requirements( requires, dependencies, store_filenames )

            mini_app = MiniAppBase.any_from_conf( coef_name, problem, cargs )
            if len(graph[coef_name]) > 1:
                for name in graph[coef_name][1:]:
                    key = 'c.' + name
                    requires.append(key)
                    dependencies[key] = getattr(coefs, name)

            problem.clear_equations()

            # Pass only the direct dependencies, not the indirect ones.
            data = {}
            for key in requires:
                data[key] = dependencies[key]

            val = mini_app(self.volume, data=data)
            setattr( coefs, coef_name, val )
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
            
        if ret_all:
            return coefs, dependencies
        else:
            return coefs

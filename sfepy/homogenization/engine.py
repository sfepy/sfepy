from sfepy.base.base import *
from sfepy.applications import SimpleApp, Application
from sfepy.fem import eval_term_op
from coefs_base import MiniAppBase

class HomogenizationEngine( SimpleApp ):

    def process_options( options ):
        get = options.get_default_attr
        
        coefs = get( 'coefs', None, 'missing "coefs" in options!' )
        requirements = get( 'requirements', None,
                            'missing "requirements" in options!' )

        return Struct( **locals() )
    process_options = staticmethod( process_options )

    def __init__( self, problem, options,
                  volume = None, output_prefix = 'he:', **kwargs ):
        """Bypasses SimpleApp.__init__()!"""
        Application.__init__( self, problem.conf, options, output_prefix,
                              **kwargs )
        self.problem = problem
        self.setup_options()
        self.setup_output_info( self.problem, self.options )
        
        if volume is None:
            self.volume = eval_term_op( None, self.app_options.total_volume,
                                        self.problem )
        else:
            self.volume = volume

    def setup_options( self ):
        SimpleApp.setup_options( self )
        po = HomogenizationEngine.process_options
        self.app_options += po( self.conf.options )
    
    def call( self, ret_all = False ):
        problem = self.problem

        opts = self.app_options
        coef_info = getattr( self.conf, opts.coefs )
        req_info = getattr( self.conf, opts.requirements )

        dependencies = {}

        coefs = Struct()

        for coef_name, cargs in coef_info.iteritems():
            output( 'computing %s...' % coef_name )
            requires = cargs.get( 'requires', [] )
            
            for req in requires:
                if req in dependencies and (dependencies[req] is not None):
                    continue

                output( 'computing dependency %s...' % req )

                rargs = req_info[req]

                save_name = rargs.get( 'save_name', '' )
                name = os.path.join( problem.output_dir, save_name )

                mini_app = MiniAppBase.any_from_conf( req, problem, rargs )
                save_hook = mini_app.make_save_hook( name,
                                                     problem.output_format,
                                                     self.post_process_hook )

#                print mini_app
                dep = mini_app( data = dependencies, save_hook = save_hook )

                dependencies[req] = dep
                output( '...done' )

#                print dep
#                pause()

            mini_app = MiniAppBase.any_from_conf( coef_name, problem, cargs )
            val = mini_app( self.volume, data = dependencies )
            print val
            setattr( coefs, coef_name, val )
#            pause()
            output( '...done' )

        if ret_all:
            return coefs, dependencies
        else:
            return coefs

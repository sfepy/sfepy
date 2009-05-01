from sfepy.base.base import *
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.applications import SimpleApp

def pde_solve(conf_filename, options=None):
    required, other = get_standard_keywords()
    conf = ProblemConf.from_file(conf_filename, required, other)
    opts = conf.options

    output_prefix = opts.get_default_attr('output_prefix', None)
    if output_prefix is None:
        output_prefix = output.prefix 

    if options is None:
        options = Struct(output_filename_trunk = None,
                         save_ebc = False,
                         save_regions = False,
                         save_field_meshes = False,
                         save_region_field_meshes = False,
                         solve_not = False)
        
    app = SimpleApp(conf, options, output_prefix)
    if hasattr( opts, 'parametric_hook' ): # Parametric study.
        parametric_hook = getattr(conf, opts.parametric_hook)
        app.parametrize(parametric_hook)

    return app()

"""
Various parametric hooks.
"""
import os
import numpy as nm

from sfepy.base.base import pause, output, default_printer
from sfepy.base.log import Log

def vary_incident_wave_dir( problem ):
    default_printer.prefix = 'vary_incident_wave_dir:'

    log_conf = {
        'is_plot' : True,
        'yscales' : ['linear'],
        'xaxes' : ['incident wave angle [degrees]'],
        'yaxes' : ['phase velocity [m/s]'],
    }
    dim = problem.domain.mesh.dim
    log = Log.from_conf( log_conf,
                         ['phase velocity %d' % ii for ii in range( dim )] )
    
    alphas = nm.linspace( 0.0, 360.0, 37 )
    output( 'running for angles:', alphas )
    pause()
    for ii, alpha in enumerate( alphas ):
        output( 'iteration %d: alpha %2f' % (ii, alpha) )
        opts = problem.conf.options

        ra = nm.pi * alpha / 180.0
        iwd = nm.zeros( (dim,), dtype = nm.float64 )
        iwd[0:2] = [nm.cos( ra ), nm.sin( ra )]
        opts.incident_wave_dir = iwd

        out = []
        yield problem, out
        
        #You have: sqrt( 10^10Pa / 10^4kg * m^3 )
        #You want:
        #Definition: 1000 m / s
        convert_to_ms = 1000.0
        phase_velocity = out[0] * convert_to_ms # to m/s.
        log( x = [alpha], *phase_velocity )

        output( 'phase velocity [m/s]:', phase_velocity )
        if opts.homogeneous:
            mat = problem.materials['matrix']
            # dilatation wave: vp = \sqrt{ (\lambda + 2\mu) / \rho }
            vd = nm.sqrt( (mat.D[0,0]) / mat.density) * convert_to_ms
            # shear wave: vs = \sqrt{ \mu / \rho }
            vs = nm.sqrt( mat.D[-1,-1] / mat.density ) * convert_to_ms
            output( 'analytical dilatation: %f, shear: %f' % (vd, vs) )
#            pause()

        yield None

    print log
    pause()

    if opts.homogeneous:
        name = 'homogeneous'
    else:
        name = 'perforated'
    fig_name = os.path.join( problem.output_dir,
                             'phase_velocity_%s%s' % (name, opts.fig_suffix) )
    log( save_figure = fig_name, finished = True )

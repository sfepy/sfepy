from sfepy.base.plotutils import pylab

from sfepy.base.base import *
from sfepy.base.la import eig
from sfepy.fem.evaluate import eval_term_op
from sfepy.base.progressbar import MyBar

##
# 26.09.2007, c
# 27.09.2007
# 01.10.2007
# 02.10.2007
def process_options( options, n_eigs ):
    try:
        save = options.save_eig_vectors
    except:
        save = (n_eigs, n_eigs)

    try:
        eig_range = options.eig_range
        if eig_range[-1] < 0:
            eig_range[-1] += n_eigs + 1
    except:
        eig_range = (0, n_eigs)
    assert eig_range[0] < (eig_range[1] - 1)
    assert eig_range[1] <= n_eigs
    
    try:
        freq_margins = 0.01 * nm.array( options.freq_margins, dtype = nm.float64 )
    except:
        freq_margins = nm.array( (0.05, 0.05), dtype = nm.float64 )

    try:
        freq_step = 0.01 * options.freq_step
    except:
        freq_step = 0.05

    try:
        feps = options.feps
    except:
        feps = 1e-8

    try:
        zeps = options.zeps
    except:
        zeps = 1e-8

    try:
        teps = options.teps
    except:
        teps = 1e-4

    try:
        eig_vector_transform = options.eig_vector_transform
    except:
        eig_vector_transform = None

    try:
        plot_tranform = options.plot_tranform
    except:
        plot_tranform = None

    try:
        squared = options.squared
    except:
        squared = True

    return Struct( **locals() )

##
# c: 08.04.2008, r: 08.04.2008
def get_method( options ):
    if hasattr( options, 'method' ):
        method = options.method
    else:
        method = 'eig.sgscipy'
    return method

##
# created:       27.09.2007
# last revision: 08.04.2008
def compute_average_density( pb, mat_name1, mat_name2 ):

    mat1 = pb.materials[mat_name1]
    region_name = mat1.region.name
    vol1 = eval_term_op( None, 'd_volume.i1.%s( u1 )' % region_name, pb )
    mat2 = pb.materials[mat_name2]
    region_name = mat2.region.name
    vol2 = eval_term_op( None, 'd_volume.i1.%s( u )' % region_name, pb )
    output( 'volumes:', vol1, vol2, vol1 + vol2 )
    output( 'densities:', mat1.density, mat2.density )

    return mat1.density * vol1 + mat2.density * vol2

##
# c: 26.09.2007, r: 08.04.2008
def compute_mass_components( pb, mtx_phi, threshold,
                           transform = None, pbar = None ):
    """Eigenmomenta..."""
    dim = pb.domain.mesh.dim
    n_dof, n_eigs = mtx_phi.shape
    n_nod = n_dof / dim
    
    term = pb.equations[0].terms[0]
    mat_name = term.get_material_names()[0]
    mat = pb.materials[mat_name]
    
    uc_name = 'uc'
    d_name = 'd'
    dot_term = 'd_volume_dot.i1.%s( %s, %s )' % (mat.region.name, d_name, uc_name)

    density = nm.empty( (n_dof,), dtype = nm.float64 )
    density.fill( mat.density )
    pb.variables[d_name].data_from_data( density, slice( 0, n_dof ) )
    
    masses = nm.empty( (n_eigs, dim), dtype = nm.float64 )

    if pbar is not None:
        pbar.init( n_eigs - 1 )
        
    for ii in xrange( n_eigs ):
        if pbar is not None:
            pbar.update( ii )
        else:
            if (ii % 100) == 0:
                output( '%d of %d (%f%%)' % (ii, n_eigs,
                                             100. * ii / (n_eigs - 1)) )
            
        if transform is None:
            vec_phi, is_zero = mtx_phi[:,ii], False
        else:
            vec_phi, is_zero = transform( mtx_phi[:,ii], (n_nod, dim) )
           
        if is_zero:
            masses[ii,:] = 0.0
        else:
            for ir in range( dim ):
                vec = vec_phi[ir::dim].copy()
                pb.variables[uc_name].data_from_data( vec, slice( 0, n_nod ) )
                val = eval_term_op( None, dot_term, pb )
                if abs( val ) >= threshold:
                    masses[ii,ir] = val
                else:
                    masses[ii,ir] = 0.0
    #            print ii, ir, val
        
    return masses

##
# 26.09.2007, c
# 27.09.2007
def compute_generalized_mass( freq, masses, eigs, average_density, squared ):
    """Assumes squared freq!"""
    dim = masses.shape[1]
    mtx_mass = nm.eye( dim, dim, dtype = nm.float64 ) * average_density
    for ir in range( dim ):
        for ic in range( dim ):
            if ir <= ic:
                if squared:
                    val = nm.sum( masses[:,ir] * masses[:,ic]\
                                  / (freq - eigs) )
                    mtx_mass[ir,ic] += - freq * val
                else:
                    val = nm.sum( masses[:,ir] * masses[:,ic]\
                                  / ((freq**2) - (eigs**2)) )
                    mtx_mass[ir,ic] += - (freq**2) * val
            else:
                mtx_mass[ir,ic] = mtx_mass[ic,ir]
    return mtx_mass


##
# c: 27.09.2007, r: 08.04.2008
def find_zero( f0, f1, masses, eigs, average_density, opts, mode ):
    feps, zeps = opts.feps, opts.zeps
    method = get_method( opts )

    fm, fp = f0, f1
    ieig = {0 : 0, 1 : -1}[mode]
    while 1:
        f = 0.5 * (fm + fp)
        mtx_mass = compute_generalized_mass( f, masses, eigs,
                                          average_density, opts.squared )
        meigs = eig( mtx_mass, eigenvectors = False, method = method )
#        print meigs

        val = meigs[ieig]
#        print f, val, fp - fm

        if (abs( val ) < zeps)\
               or ((fp - fm) < (100.0 * nm.finfo( float ).eps))\
               or ((fp - fm) < feps):
            return 0, f, val

        if mode == 0:
            if (f - f0) < feps:
                return 2, f0, val
            elif (f1 - f) < feps:
                return 1, f, val
        elif mode == 1:
            if (f1 - f) < feps:
                return 1, f, val
            elif (f - f0) < feps:
                return 2, f0, val
            
        if val > 0.0:
            fp = f
        else:
            fm = f

##
# c: 27.09.2007, r: 08.04.2008
def describe_gaps( gaps ):
    kinds = []
    for ii, (gmin, gmax) in enumerate( gaps ):

        if (gmin[0] == 2) and (gmax[0] == 2):
            kind = ('p', 'propagation zone')
        elif (gmin[0] == 1) and (gmax[0] == 2):
            kind = ('w', 'full weak band gap')
        elif (gmin[0] == 0) and (gmax[0] == 2):
            kind = ('wp', 'weak band gap + propagation zone')
        elif (gmin[0] == 1) and (gmax[0] == 1):
            kind = ('s', 'full strong band gap (due to end of freq. range or'
                    ' too large thresholds)')
        elif (gmin[0] == 1) and (gmax[0] == 0):
            kind = ('sw', 'strong band gap + weak band gap')
        elif (gmin[0] == 0) and (gmax[0] == 0):
            kind = ('swp', 'strong band gap + weak band gap + propagation zone')
        else:
            output( 'impossible band gap combination:' )
            output( gmin, gmax )
            raise ValueError
        kinds.append( kind )
    return kinds

##
# created:       01.10.2007
# last revision: 13.12.2007
def transform_plot_data( datas, plot_tranform, funmod ):
    if plot_tranform is not None:
        fun = getattr( funmod, plot_tranform[0] )

    dmin, dmax = 1e+10, -1e+10
    tdatas = []
    for data in datas:
        tdata = data.copy()
        if plot_tranform is not None:
            tdata[:,1:] = fun( tdata[:,1:], *plot_tranform[1:] )
        dmin = min( dmin, tdata[:,1:].min() )
        dmax = max( dmax, tdata[:,1:].max() )
        tdatas.append( tdata )
    dmin, dmax = min( dmax - 1e-8, dmin ), max( dmin + 1e-8, dmax )
    return (dmin, dmax), tdatas

##
# c: 27.09.2007, r: 12.06.2008
def plot_logs( fig_num, logs, freq_range, plot_range, squared, show = False ):
    if pylab is None: return

    fig = pylab.figure( fig_num )
    ax = fig.add_subplot( 111 )

    for f in freq_range:
        l0 = ax.plot( [f, f], plot_range, 'r' )

    for log in logs:
        l1 = ax.plot( log[:,0], log[:,1], 'b--' )
        l2 = ax.plot( log[:,0], log[:,2], 'b-' )

    fmin, fmax = logs[0][0,0], logs[-1][-1,0]
    ax.plot( [fmin, fmax], [0, 0], 'k--' )
    ax.legend( (l0, l1, l2),
               ('eigenfrequencies', 'min eig($a^*$)', 'max eig($a^*$)') )
    if squared:
        ax.set_xlabel( r'$\lambda$, $\omega^2$' )
    else:
        ax.set_xlabel( r'$\sqrt{\lambda}$, $\omega$' )
    ax.set_ylabel( r'eigenvalues of mass matrix $A^*$' )

    if show:
        ax.set_xlim( [fmin, fmax] )
        ax.set_ylim( plot_range )
        pylab.show()
    
##
# c: 27.09.2007, r: 12.06.2008
def plot_gaps( fig_num, gaps, kinds, freq_range, plot_range, show = False ):
    if pylab is None: return

    def draw_rect( ax, x, y, color ):
        ax.fill( nm.asarray( x )[[0,1,1,0]],
                 nm.asarray( y )[[0,0,1,1]],
                 fc = color, linewidth = 0 )

    fig = pylab.figure( fig_num )
    ax = fig.add_subplot( 111 )

    # Colors.
    strong = (1, 1, 0.5)
    weak = (1, 1, 1)
    propagation = (0.5, 1, 0.5)

    for ii in xrange( len( freq_range ) - 1 ):
        f0, f1 = freq_range[[ii, ii+1]]
        gmin, gmax = gaps[ii]
        kind, kind_desc = kinds[ii]

        if kind == 'p':
            draw_rect( ax, (f0, f1), plot_range, propagation )
            info = [(f0, f1)]
        elif kind == 'w':
            draw_rect( ax, (f0, f1), plot_range, weak )
            info = [(f0, f1)]
        elif kind == 'wp':
            draw_rect( ax, (f0, gmin[1]), plot_range, weak )
            draw_rect( ax, (gmin[1], f1), plot_range, propagation )
            info = [(f0, gmin[1]), (gmin[1], f1)]
        elif kind == 's':
            draw_rect( ax, (f0, f1), plot_range, strong )
            info = [(f0, f1)]
        elif kind == 'sw':
            draw_rect( ax, (f0, gmax[1]), plot_range, strong )
            draw_rect( ax, (gmax[1], f1), plot_range, weak )
            info = [(f0, gmax[1]), (gmax[1], f1)]
        elif kind == 'swp':
            draw_rect( ax, (f0, gmax[1]), plot_range, strong )
            draw_rect( ax, (gmax[1], gmin[1]), plot_range, weak )
            draw_rect( ax, (gmin[1], f1), plot_range, propagation )
            info = [(f0, gmax[1]), (gmax[1], gmin[1]), (gmin[1], f1)]
        else:
            output( 'impossible band gap combination:' )
            output( gmin, gmax )
            raise ValueError

        output( ii, gmin[0], gmax[0], '%.8f' % f0, '%.8f' % f1 )
        output( ' -> %s\n    %s' %(kind_desc, info) )

    if show:
        ax.set_xlim( [freq_range[0], freq_range[-1]] )
        ax.set_ylim( plot_range )
        pylab.show()
    
##
# c: 27.09.2007, r: 08.04.2008
def detect_band_gaps( pb, eigs, mtx_phi, conf, options ):
    
    average_density = compute_average_density( pb, 'matrix', 'inclusion' )
    output( 'average density:', average_density )

    n_eigs = eigs.shape[0]
    opts = process_options( conf.options, n_eigs )
    method = get_method( conf.options )
    output( 'method:', method )
    
    if not opts.squared:
        eigs = nm.sqrt( eigs )

    freq_range = eigs[slice( *opts.eig_range )]
    n_freq = freq_range.shape[0]
    min_freq, max_freq = freq_range[0], freq_range[-1]
    margins = opts.freq_margins * (max_freq - min_freq)

    prev_eig = min_freq - margins[0]
    next_eig = max_freq + margins[1]
    if opts.eig_range[0] > 0:
        prev_eig = max( eigs[opts.eig_range[0]-1] + opts.feps, prev_eig )
    if opts.eig_range[1] < n_eigs:
        next_eig = min( eigs[opts.eig_range[1]] - opts.feps, next_eig )
    prev_eig = max( opts.feps, prev_eig )
    next_eig = max( opts.feps, next_eig, prev_eig + opts.feps )
    freq_range_margins = nm.r_[prev_eig, freq_range, next_eig]

    output( 'freq. range             : [%8.3f, %8.3f]' % (min_freq, max_freq) )
    output( 'freq. range with margins: [%8.3f, %8.3f]'\
          % tuple( freq_range_margins[[0,-1]] ) )

##     print freq_range
##     print freq_range_margins
##     pause()

    if opts.eig_vector_transform is not None:
        fun = getattr( conf.funmod, opts.eig_vector_transform[0] )
        def _wrap_transform( vec, shape ):
            return fun( vec, shape, *opts.eig_vector_transform[1:] )
    else:
        _wrap_transform = None
    output( 'mass matrix components...')
    pbar = MyBar( 'computing:' )
    tt = time.clock()
    masses = compute_mass_components( pb, mtx_phi, opts.teps, _wrap_transform, pbar )
    output( '...done in %.2f s' % (time.clock() - tt) )

    logs = []
    gaps = []
    df = opts.freq_step * (max_freq - min_freq)
    for ii in xrange( n_freq + 1 ):
        f0, f1 = freq_range_margins[[ii, ii+1]]
        output( 'interval: ]%.8f, %.8f[...' % (f0, f1) )

        log = []
        num = max( 5, (f1 - f0) / df )
        log_freqs = nm.linspace( f0 + opts.feps, f1 - opts.feps, num )
        for f in log_freqs:
            mtx_mass = compute_generalized_mass( f, masses, eigs,
                                              average_density, opts.squared )
            meigs = eig( mtx_mass, eigenvectors = False, method = method )
            log.append( [f, meigs[0], meigs[-1]] )
        log0, log1 = log[0], log[-1]
        if log0[1] > 0.0: # No gaps.
            gap = ([2, f0, log0[1]], [2, f0, log0[2]])
        elif log1[2] < 0.0: # Full interval strog gap.
            gap = ([1, f1, log1[1]], [1, f1, log1[2]])
        else:
            output( 'finding zero of the largest eig...' )
            smax, fmax, vmax = find_zero( f0, f1, masses, eigs,
                                         average_density, opts, 1 )
            output( '...done' )
            if smax in [0, 2]:
                output( 'finding zero of the smallest eig...' )
                smin, fmin, vmin = find_zero( fmax, f1, masses, eigs,
                                             average_density, opts, 0 )
                output( '...done' )
            elif smax == 1:
                smin = 1 # both are negative everywhere.
                fmin, vmin = fmax, vmax

            gap = ([smin, fmin, vmin], [smax, fmax, vmax])

        output( gap[0] )
        output( gap[1] )
#        pause()
        gaps.append( gap )
        logs.append( nm.array( log, dtype = nm.float64 ) )
        output( '...done' )

    kinds = describe_gaps( gaps )

    return Struct( logs = logs, gaps = gaps, kinds = kinds,
                   freq_range = freq_range, freq_range_margins = freq_range_margins,
                   opts = opts )

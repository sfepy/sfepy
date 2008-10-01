from sfepy.base.plotutils import pylab

from sfepy.base.base import *
from sfepy.base.la import eig
from sfepy.fem.evaluate import eval_term_op
from sfepy.base.progressbar import MyBar

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
        fixed_eig_range = options.fixed_eig_range
    except:
        fixed_eig_range = None

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
        teps_rel = options.teps_rel
    except:
        teps_rel = True

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

    try:
        plot_options = options.plot_options
    except:
        plot_options = {
            'show' : True,
            'legend' : False,
        }

    try:
        fig_name = options.fig_name
    except:
        fig_name = None

    default_plot_rsc =  {
        'resonance' : {'linewidth' : 0.5, 'color' : 'r', 'linestyle' : '-' },
        'masked' : {'linewidth' : 0.5, 'color' : 'r', 'linestyle' : ':' },
        'x_axis' : {'linewidth' : 0.5, 'color' : 'k', 'linestyle' : '--' },
        'eig_min' : {'linewidth' : 0.5, 'color' : 'b', 'linestyle' : '--' },
        'eig_max' : {'linewidth' : 0.5, 'color' : 'b', 'linestyle' : '-' },
        'strong_gap' : {'linewidth' : 0, 'facecolor' : (1, 1, 0.5) },
        'weak_gap' : {'linewidth' : 0, 'facecolor' : (1, 1, 1) },
        'propagation' : {'linewidth' : 0, 'facecolor' : (0.5, 1, 0.5) },
        'params' : {'axes.labelsize': 'large',
                    'text.fontsize': 'large',
                    'legend.fontsize': 'large',
                    'xtick.labelsize': 'large',
                    'ytick.labelsize': 'large',
                    'text.usetex': False},
    }
    try:
        plot_rsc = options.plot_rsc
        # Missing values are set to default.
        for key, val in default_plot_rsc.iteritems():
            if not key in plot_rsc:
                plot_rsc[key] = val
    except:
        plot_rsc = default_plot_rsc
    del default_plot_rsc

    try:
        eigenmomentum = options.eigenmomentum
    except:
        raise ValueError( 'missing key "eigenmomentum" in options!' )

    try:
        region_to_material = options.region_to_material
    except:
        raise ValueError( 'missing key "region_to_material" in options!' )

    try:
        volume = options.volume
    except:
        raise ValueError( 'missing key "volume" in options!' )
    
    return Struct( **locals() )

##
# c: 08.04.2008, r: 08.04.2008
def get_method( options ):
    if hasattr( options, 'method' ):
        method = options.method
    else:
        method = 'eig.sgscipy'
    return method

def compute_average_density( pb, volume_term, region_to_material ):

    average_density = 0.0
    total_volume = 0.0
    for region_name, mat_name in region_to_material.iteritems():
        mat = pb.materials[mat_name]
        assert region_name == mat.region.name
        vol = eval_term_op( None, volume_term % region_name, pb )
        density = mat.get_data( region_name, mat.igs[0], 'density' )
        output( 'region %s: volume %f, density %f' % (region_name,
                                                      vol, density ) )
        average_density += vol * density
        total_volume += vol
    output( 'total volume:', total_volume )

    average_density /= total_volume

    return average_density

def compute_eigenmomenta( pb, conf_eigenmomentum, region_to_material,
                          eig_vectors, threshold, threshold_is_relative,
                          transform = None, pbar = None ):
    """Eigenmomenta.

    Valid == True means an eigenmomentum above threshold."""
    dim = pb.domain.mesh.dim
    n_dof, n_eigs = eig_vectors.shape
    n_nod = n_dof / dim

    u_name = conf_eigenmomentum['var']
    rnames = conf_eigenmomentum['regions']
    term = conf_eigenmomentum['term']
    tt = []
    for rname in rnames:
        mat = pb.materials[region_to_material[rname]]
        density = mat.get_data( rname, mat.igs[0], 'density' )
        tt.append( term % (density, rname, u_name) )
    em_eq = ' + '.join( tt )
    output( 'eigenmomentum equation:', em_eq )

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
            vec_phi, is_zero = eig_vectors[:,ii], False
        else:
            vec_phi, is_zero = transform( eig_vectors[:,ii], (n_nod, dim) )
           
        if is_zero:
            masses[ii,:] = 0.0
        else:
            pb.variables[u_name].data_from_data( vec_phi.copy() )
            val = eval_term_op( None, em_eq, pb )
            masses[ii,:] = val
#            print ii, ir, val

    mag = nm.zeros( (n_eigs,), dtype = nm.float64 )
    for ir in range( dim ):
        mag += masses[:,ir] ** 2
    mag = nm.sqrt( mag )

    if threshold_is_relative:
        tol = threshold * mag.max()
    else:
        tol = threshold
        
    valid = nm.where( mag < tol, False, True )
    mask = nm.where( valid == False )[0]
    masses[mask,:] = 0.0
    n_zeroed = mask.shape[0]

    output( '%d of %d eigenmomenta zeroed (under %.2e)'\
            % (n_zeroed, n_eigs, tol) )
##     print valid
##     import pylab
##     pylab.plot( masses )
##     pylab.show()
    
    return n_zeroed, valid, masses

def compute_generalized_mass( freq, masses, eigs, average_density, squared ):
    """Assumes that `masses`, `eigs` contain only valid resonances."""
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

def find_zero( f0, f1, masses, eigs, average_density, opts, mode ):
    """For f \in ]f0, f1[ find frequency f for which either the smallest
    (`mode` = 0) or the largest (`mode` = 1) eigenvalue of M is zero.

    Return:

    (flag, frequency, eigenvalue)

    mode | flag | meaning
    0, 1 | 0    | eigenvalue -> 0 for f \in ]f0, f1[
    0    | 1    | f -> f1, smallest eigenvalue < 0
    0    | 2    | f -> f0, smallest eigenvalue > 0 and -> -\infty
    1    | 1    | f -> f1, largest eigenvalue < 0 and  -> +\infty
    1    | 2    | f -> f0, largest eigenvalue > 0
    """
    feps, zeps = opts.feps, opts.zeps
    method = get_method( opts )

    fm, fp = f0, f1
    ieig = {0 : 0, 1 : -1}[mode]
    while 1:
        f = 0.5 * (fm + fp)
        mtx_mass = compute_generalized_mass( f, masses, eigs,
                                             average_density, opts.squared )
        meigs = eig( mtx_mass, eigenvectors = False, method = method )
##         print meigs

        val = meigs[ieig]
##         print f, f0, f1, fm, fp, val

        if (abs( val ) < zeps)\
               or ((fp - fm) < (100.0 * nm.finfo( float ).eps)):
            return 0, f, val

        if mode == 0:
            if (f - f0) < feps:
                return 2, f0, val
            elif (f1 - f) < feps:
                return 1, f1, val
        elif mode == 1:
            if (f1 - f) < feps:
                return 1, f1, val
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

def plot_eigs( fig_num, plot_rsc, valid, freq_range, plot_range,
               show = False, clear = False, new_axes = False ):
    """
    Plot resonance/eigen-frequencies.

    `valid` must correspond to `freq_range`

    resonances : red
    masked resonances: dotted red
    """
    if pylab is None: return
    assert len( valid ) == len( freq_range )

    fig = pylab.figure( fig_num )
    if clear:
        fig.clf()
    if new_axes:
        ax = fig.add_subplot( 111 )
    else:
        ax = fig.gca()

    l0 = l1 = None
    for ii, f in enumerate( freq_range ):
        if valid[ii]:
            l0 = ax.plot( [f, f], plot_range, **plot_rsc['resonance'] )[0]
        else:
            l1 = ax.plot( [f, f], plot_range, **plot_rsc['masked'] )[0]
 
    if l0:
        l0.set_label( 'eigenfrequencies' )
    if l1:
        l1.set_label( 'masked eigenfrequencies' )

    if new_axes:
        ax.set_xlim( [freq_range[0], freq_range[-1]] )
        ax.set_ylim( plot_range )

    if show:
        pylab.show()
    return fig 

def plot_logs( fig_num, plot_rsc,
               logs, valid, freq_range, plot_range, squared,
               draw_eigs = True, show_legend = True, show = False,
               clear = False, new_axes = False ):
    """
    Plot logs of min/max eigs of M.
    """
    if pylab is None: return

    fig = pylab.figure( fig_num )
    if clear:
        fig.clf()
    if new_axes:
        ax = fig.add_subplot( 111 )
    else:
        ax = fig.gca()

    if draw_eigs:
        aux = plot_eigs( fig_num, plot_rsc, valid, freq_range, plot_range )

    for log in logs:
        l1 = ax.plot( log[:,0], log[:,1], **plot_rsc['eig_min'] )
        l2 = ax.plot( log[:,0], log[:,2], **plot_rsc['eig_max'] )
    l1[0].set_label( 'min eig($M^*$)' )
    l2[0].set_label( 'max eig($M^*$)' )

    fmin, fmax = logs[0][0,0], logs[-1][-1,0]
    ax.plot( [fmin, fmax], [0, 0], **plot_rsc['x_axis'] )

    if squared:
        ax.set_xlabel( r'$\lambda$, $\omega^2$' )
    else:
        ax.set_xlabel( r'$\sqrt{\lambda}$, $\omega$' )
    ax.set_ylabel( r'eigenvalues of mass matrix $M^*$' )

    if new_axes:
        ax.set_xlim( [fmin, fmax] )
        ax.set_ylim( plot_range )

    if show_legend:
        ax.legend()

    if show:
        pylab.show()
    return fig
    
def plot_gaps( fig_num, plot_rsc, gaps, kinds, freq_range, plot_range,
               show = False, clear = False, new_axes = False ):
    """ """
    if pylab is None: return

    def draw_rect( ax, x, y, rsc ):
        ax.fill( nm.asarray( x )[[0,1,1,0]],
                 nm.asarray( y )[[0,0,1,1]],
                 **rsc )

    fig = pylab.figure( fig_num )
    if clear:
        fig.clf()
    if new_axes:
        ax = fig.add_subplot( 111 )
    else:
        ax = fig.gca()

    # Colors.
    strong = plot_rsc['strong_gap']
    weak = plot_rsc['weak_gap']
    propagation = plot_rsc['propagation']

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

    if new_axes:
        ax.set_xlim( [freq_range[0], freq_range[-1]] )
        ax.set_ylim( plot_range )

    if show:
        pylab.show()
    return fig

def cut_freq_range( freq_range, eigs, valid, freq_margins, eig_range,
                    fixed_eig_range, feps ):
    """Cut off masked resonance frequencies. Margins are preserved, like no
    resonances were cut.

    Return:

      freq_range - new resonance frequencies
      freq_range_margins - freq_range with prepended/appended margins
                           equal to fixed_eig_range if it is not None
    """
    n_freq = freq_range.shape[0]
    n_eigs = eigs.shape[0]

    output( 'masked resonance frequencies in range:' )
    valid_slice = slice( *eig_range )
    output( nm.where( valid[valid_slice] == False )[0] )

    if fixed_eig_range is None:
        min_freq, max_freq = freq_range[0], freq_range[-1]
        margins = freq_margins * (max_freq - min_freq)
        prev_eig = min_freq - margins[0]
        next_eig = max_freq + margins[1]
        if eig_range[0] > 0:
            prev_eig = max( eigs[eig_range[0]-1] + feps, prev_eig )
        if eig_range[1] < n_eigs:
            next_eig = min( eigs[eig_range[1]] - feps, next_eig )
        prev_eig = max( feps, prev_eig )
        next_eig = max( feps, next_eig, prev_eig + feps )
    else:
        prev_eig, next_eig = fixed_eig_range

    freq_range = freq_range[valid[valid_slice]]
    freq_range_margins = nm.r_[prev_eig, freq_range, next_eig]

    return freq_range, freq_range_margins
    
def detect_band_gaps( pb, eigs, eig_vectors, conf, options ):
    """Detect band gaps given solution to eigenproblem (eigs,
    eig_vectors). Only valid resonance frequencies (e.i. those for which
    corresponding eigenmomenta are above a given threshold) are taken into
    account.

    ??? make feps relative to ]f0, f1[ size ???
    """
    n_eigs = eigs.shape[0]
    opts = process_options( conf.options, n_eigs )
    method = get_method( conf.options )
    output( 'method:', method )

    average_density = compute_average_density( pb,
                                               opts.volume,
                                               opts.region_to_material )
    output( 'average density:', average_density )
    
    if not opts.squared:
        eigs = nm.sqrt( eigs )

    if opts.fixed_eig_range is not None:
        mine, maxe = opts.fixed_eig_range
        ii = nm.where( (eigs > mine) & (eigs < maxe) )[0]
        freq_range_initial = eigs[ii]
        opts.eig_range = (ii[0], ii[-1]+1) # +1 as it is a slice.
    else:
        freq_range_initial = eigs[slice( *opts.eig_range )]
    output( 'initial freq. range     : [%8.3f, %8.3f]'\
            % tuple( freq_range_initial[[0,-1]] ) )

    if opts.eig_vector_transform is not None:
        fun = getattr( conf.funmod, opts.eig_vector_transform[0] )
        def wrap_transform( vec, shape ):
            return fun( vec, shape, *opts.eig_vector_transform[1:] )
    else:
        wrap_transform = None
    output( 'mass matrix eigenmomenta...')
    pbar = MyBar( 'computing:' )
    tt = time.clock()
    aux = compute_eigenmomenta( pb, opts.eigenmomentum,
                                opts.region_to_material,
                                eig_vectors, opts.teps, opts.teps_rel,
                                wrap_transform, pbar )
    n_zeroed, valid, masses = aux
    output( '...done in %.2f s' % (time.clock() - tt) )
    aux = cut_freq_range( freq_range_initial, eigs, valid,
                          opts.freq_margins, opts.eig_range,
                          opts.fixed_eig_range,
                          opts.feps )
    freq_range, freq_range_margins = aux
    if len( freq_range ):
        output( 'freq. range             : [%8.3f, %8.3f]'\
                % tuple( freq_range[[0,-1]] ) )
    else:
        # All masked.
        output( 'freq. range             : all masked!' )

    min_freq, max_freq = freq_range_margins[0], freq_range_margins[-1]
    output( 'freq. range with margins: [%8.3f, %8.3f]'\
            % (min_freq, max_freq) )
        
    logs = []
    gaps = []

    df = opts.freq_step * (max_freq - min_freq)
    valid_masses = masses[valid,:]
    valid_eigs = eigs[valid]
    cgm = compute_generalized_mass
    for ii in xrange( freq_range.shape[0] + 1 ):

        f0, f1 = freq_range_margins[[ii, ii+1]]
        output( 'interval: ]%.8f, %.8f[...' % (f0, f1) )

        log = []
        num = min( 200, max( 20, (f1 - f0) / df ) )
        log_freqs = nm.linspace( f0 + opts.feps, f1 - opts.feps, num )
        for f in log_freqs:
            mtx_mass = cgm( f, valid_masses, valid_eigs,
                            average_density, opts.squared )
            meigs = eig( mtx_mass, eigenvectors = False, method = method )
            log.append( [f, meigs[0], meigs[-1]] )

        log0, log1 = log[0], log[-1]
        if log0[1] > 0.0: # No gaps.
            gap = ([2, f0, log0[1]], [2, f0, log0[2]])
        elif log1[2] < 0.0: # Full interval strong gap.
            gap = ([1, f1, log1[1]], [1, f1, log1[2]])
        else:
            # Insert fmin, fmax into log.
            alog = nm.array( log )

            output( 'finding zero of the largest eig...' )
            smax, fmax, vmax = find_zero( f0, f1, valid_masses, valid_eigs,
                                          average_density, opts, 1 )
            mtx_mass = cgm( fmax, valid_masses, valid_eigs,
                            average_density, opts.squared )
            meigs = eig( mtx_mass, eigenvectors = False, method = method )
            im = nm.searchsorted( alog[:,0], fmax )
            log.insert( im, (fmax, meigs[0], meigs[-1] ) )

            output( '...done' )
            if smax in [0, 2]:
                output( 'finding zero of the smallest eig...' )
                # having fmax instead of f0 does not work if feps is large.
                smin, fmin, vmin = find_zero( f0, f1, valid_masses, valid_eigs,
                                              average_density, opts, 0 )
                mtx_mass = cgm( fmin, valid_masses, valid_eigs,
                                average_density, opts.squared )
                meigs = eig( mtx_mass, eigenvectors = False, method = method )
                im = nm.searchsorted( alog[:,0], fmin )
                # +1 due to fmax already inserted before.
                log.insert( im+1, (fmin, meigs[0], meigs[-1] ) )

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
                   valid = valid, eig_range = slice( *opts.eig_range ),
                   n_eigs = eigs.shape[0], n_zeroed = n_zeroed,
                   freq_range_initial = freq_range_initial,
                   freq_range = freq_range,
                   freq_range_margins = freq_range_margins,
                   opts = opts )

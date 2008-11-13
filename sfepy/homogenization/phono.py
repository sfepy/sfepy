from sfepy.base.plotutils import pylab

from sfepy.base.base import *
from sfepy.base.la import eig
from sfepy.fem.evaluate import eval_term_op
from sfepy.base.progressbar import MyBar
from sfepy.homogenization.utils import coor_to_sym

class AcousticMassTensor( Struct ):

    def __init__( self, eigenmomenta, eigs, dv_info ):
        Struct.__init__( self, eigenmomenta = eigenmomenta,
                         eigs = eigs, dv_info = dv_info )

    def __call__( self, freq ):
        """`eigenmomenta`, `eigs` should contain only valid resonances."""
        ema = self.eigenmomenta 
        eigs = self.eigs

        dim = ema.shape[1]
        fmass = nm.zeros( (dim, dim), dtype = nm.float64 )

        de = (freq**2) - (eigs)
        if not nm.isfinite( de ).all():
            raise ValueError( 'frequency %e too close to resonance!' % freq )

        for ir in range( dim ):
            for ic in range( dim ):
                if ir <= ic:
                    val = nm.sum( ema[:,ir] * ema[:,ic] / de )
                    fmass[ir,ic] += (freq**2) * val
                else:
                    fmass[ir,ic] = fmass[ic,ir]

        eye = nm.eye( dim, dim, dtype = nm.float64 )
        mtx_mass = (eye * self.dv_info.average_density) \
                   - (fmass / self.dv_info.total_volume)

        return mtx_mass

class AppliedLoadTensor( Struct ):

    def __init__( self, eigenmomenta, ueigenmomenta, eigs, dv_info ):
        Struct.__init__( self, eigenmomenta = eigenmomenta,
                         eigs = eigs, dv_info = dv_info )


    def __call__( self, freq ):
        """`eigenmomenta`, `ueigenmomenta`, `eigs` should contain only valid
        resonances."""
        ema, uema = self.eigenmomenta, self.ueigenmomenta

        dim = ema.shape[1]
        fload = nm.zeros( (dim, dim), dtype = nm.float64 )

        de = (freq**2) - (eigs)
        if not nm.isfinite( de ).all():
            raise ValueError( 'frequency %e too close to resonance!' % freq )

        for ir in range( dim ):
            for ic in range( dim ):
                val = nm.sum( ema[:,ir] * uema[:,ic] / de )
                fload[ir,ic] += (freq**2) * val

        eye = nm.eye( (dim, dim), dtype = nm.float64 )

        mtx_load = eye - (fload / self.dv_info.total_volume)

        return mtx_load

def get_callback( mass, method, christoffel = None, mode = 'trace' ):
    """
    Return callback to solve band gaps or dispersion eigenproblem P.

    Find zero callbacks return:
      eigenvalues

    Trace callbacks return:
      (eigenvalues, )
    or 
      (eigenvalues, eigenvectors) (in full (dispoersion) mode)

    If christoffel is None, P is
      M w = \lambda w,
    otherwise it is
      omega^2 M w = \eta \Gamma w"""

    def find_zero_callback( f ):
        meigs = eig( mass( f ), eigenvectors = False, method = method )
        return meigs

    def find_zero_full_callback( f ):
        meigs = eig( (f**2) * mass( f ), mtx_b = christoffel,
                     eigenvectors = False, method = method )
        return meigs

    def trace_callback( f ):
        meigs = eig( mass( f ), eigenvectors = False, method = method )
        return meigs,

    def trace_full_callback( f ):
        meigs, mvecs = eig( (f**2) * mass( f ), mtx_b = christoffel,
                            eigenvectors = True, method = method )
        
        return meigs, mvecs

    if christoffel is not None:
        mode += '_full'

    return eval( mode + '_callback' )

def compute_density_volume_info( pb, volume_term, region_to_material ):
    """Computes volumes, densities of regions specified in
    `region_to_material`, average density and total volume."""
    average_density = 0.0
    total_volume = 0.0
    volumes = {}
    densities = {}
    for region_name, mat_name in region_to_material.iteritems():
        mat = pb.materials[mat_name]
        assert_( region_name == mat.region.name )
        vol = eval_term_op( None, volume_term % region_name, pb )
        density = mat.get_data( region_name, mat.igs[0], 'density' )
        output( 'region %s: volume %f, density %f' % (region_name,
                                                      vol, density ) )

        volumes[region_name] = vol
        densities[region_name] = density

        average_density += vol * density
        total_volume += vol
    output( 'total volume:', total_volume )

    average_density /= total_volume

    return Struct( name = 'density_volume_info',
                   average_density = average_density,
                   total_volume = total_volume,
                   volumes = volumes,
                   densities = densities )

def compute_eigenmomenta( em_equation, u_name, pb, eig_vectors,
                          transform = None, pbar = None ):

    dim = pb.domain.mesh.dim
    n_dof, n_eigs = eig_vectors.shape

    eigenmomenta = nm.empty( (n_eigs, dim), dtype = nm.float64 )

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
            eigenmomenta[ii,:] = 0.0

        else:
            pb.variables[u_name].data_from_data( vec_phi.copy() )
            val = eval_term_op( None, em_equation, pb )
            eigenmomenta[ii,:] = val

    return eigenmomenta

def prepare_eigenmomenta( pb, conf_eigenmomentum, region_to_material,
                          eig_vectors, threshold, threshold_is_relative,
                          unweighted = False, transform = None, pbar = None ):
    """Eigenmomenta.

    unweighted = True: compute also unweighted eigenmomenta needed for applied
    load tensor
    
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
    output( 'equation:', em_eq )

    eigenmomenta = compute_eigenmomenta( em_eq, u_name, pb, eig_vectors,
                                         transform, pbar )

    mag = nm.zeros( (n_eigs,), dtype = nm.float64 )
    for ir in range( dim ):
        mag += eigenmomenta[:,ir] ** 2
    mag = nm.sqrt( mag )

    if threshold_is_relative:
        tol = threshold * mag.max()
    else:
        tol = threshold
        
    valid = nm.where( mag < tol, False, True )
    mask = nm.where( valid == False )[0]
    eigenmomenta[mask,:] = 0.0
    n_zeroed = mask.shape[0]

    output( '%d of %d eigenmomenta zeroed (under %.2e)'\
            % (n_zeroed, n_eigs, tol) )
##     print valid
##     import pylab
##     pylab.plot( eigenmomenta )
##     pylab.show()

    if unweighted:
        tt = []
        for rname in rnames:
            tt.append( term % (1.0, rname, u_name) )
        uem_eq = ' + '.join( tt )
        output( 'unweighted:', uem_eq )

        ueigenmomenta = compute_eigenmomenta( uem_eq, u_name, pb, eig_vectors,
                                     transform, pbar )
        ueigenmomenta[mask,:] = 0.0
        return n_zeroed, valid, eigenmomenta, ueigenmomenta

    else:
        return n_zeroed, valid, eigenmomenta

def compute_cat( mtx_d, iw_dir ):
    """Compute Christoffel acoustic tensor given the elasticity tensor and
    incident wave direction (unit vector).

    \Gamma_{ik} = D_{ijkl} n_j n_l
    """
    dim = iw_dir.shape[0]

    cat = nm.zeros( (dim, dim), dtype = nm.float64 )
    for ii in range( dim ):
        for ij in range( dim ):
            ir = coor_to_sym( ii, ij, dim )
            for ik in range( dim ):
                for il in range( dim ):
                    ic = coor_to_sym( ik, il, dim )
                    cat[ii,ik] += mtx_d[ir,ic] * iw_dir[ij] * iw_dir[il]

    return cat

def compute_polarization_angles( iw_dir, wave_vectors ):
    """Computes angle between incident wave direction `iw_dir` and wave
    vectors. Vector length does not matter (can use eigenvectors directly)."""
    pas = []

    iw_dir = iw_dir / nla.norm( iw_dir )
    idims = range( iw_dir.shape[0] )
    for vecs in wave_vectors:
        pa = nm.empty( vecs.shape[:-1], dtype = nm.float64 )
        for ir, vec in enumerate( vecs ):
            for ic in idims:
                vv = vec[:,ic]
                pa[ir,ic] = nm.arccos( nm.dot( iw_dir, vv ) / nla.norm( vv ) )

        pas.append( pa )

    return pas

def find_zero( f0, f1, callback, feps, zeps, mode ):
    """
    For f \in ]f0, f1[ find frequency f for which either the smallest (`mode` =
    0) or the largest (`mode` = 1) eigenvalue of problem P given by `callback`
    is zero.

    Return:

    (flag, frequency, eigenvalue)

    mode | flag | meaning
    0, 1 | 0    | eigenvalue -> 0 for f \in ]f0, f1[
    0    | 1    | f -> f1, smallest eigenvalue < 0
    0    | 2    | f -> f0, smallest eigenvalue > 0 and -> -\infty
    1    | 1    | f -> f1, largest eigenvalue < 0 and  -> +\infty
    1    | 2    | f -> f0, largest eigenvalue > 0
    """
    fm, fp = f0, f1
    ieig = {0 : 0, 1 : -1}[mode]
    while 1:
        f = 0.5 * (fm + fp)
        meigs = callback( f )
##         print meigs

        val = meigs[ieig]
##         print f, f0, f1, fm, fp, val
##         print '%.16e' % f, '%.16e' % fm, '%.16e' % fp, '%.16e' % val

        if (abs( val ) < zeps)\
               or ((fp - fm) < (abs( fm ) * nm.finfo( float ).eps)):
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
            msg = 'impossible band gap combination: %d, %d' % (gmin, gmax)
            raise ValueError( msg )
        kinds.append( kind )
    return kinds

def transform_plot_data( datas, plot_transform, funmod ):
    if plot_transform is not None:
        fun = getattr( funmod, plot_transform[0] )

    dmin, dmax = 1e+10, -1e+10
    tdatas = []
    for data in datas:
        tdata = data.copy()
        if plot_transform is not None:
            tdata = fun( tdata, *plot_transform[1:] )
        dmin = min( dmin, tdata.min() )
        dmax = max( dmax, tdata.max() )
        tdatas.append( tdata )
    dmin, dmax = min( dmax - 1e-8, dmin ), max( dmin + 1e-8, dmax )
    return (dmin, dmax), tdatas

def plot_eigs( fig_num, plot_rsc, plot_labels, valid, freq_range, plot_range,
               show = False, clear = False, new_axes = False ):
    """
    Plot resonance/eigen-frequencies.

    `valid` must correspond to `freq_range`

    resonances : red
    masked resonances: dotted red
    """
    if pylab is None: return
    assert_( len( valid ) == len( freq_range ) )

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
        l0.set_label( plot_labels['resonance'] )
    if l1:
        l1.set_label( plot_labels['masked'] )

    if new_axes:
        ax.set_xlim( [freq_range[0], freq_range[-1]] )
        ax.set_ylim( plot_range )

    if show:
        pylab.show()
    return fig 

def plot_logs( fig_num, plot_rsc, plot_labels,
               freqs, logs, valid, freq_range, plot_range, squared,
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
        aux = plot_eigs( fig_num, plot_rsc, plot_labels, valid, freq_range,
                         plot_range )

    for ii, log in enumerate( logs ):
        l1 = ax.plot( freqs[ii], log[:,0], **plot_rsc['eig_min'] )
        l2 = ax.plot( freqs[ii], log[:,-1], **plot_rsc['eig_max'] )
    l1[0].set_label( plot_labels['eig_min'] )
    l2[0].set_label( plot_labels['eig_max'] )

    fmin, fmax = freqs[0][0], freqs[-1][-1]
    ax.plot( [fmin, fmax], [0, 0], **plot_rsc['x_axis'] )

    if squared:
        ax.set_xlabel( r'$\lambda$, $\omega^2$' )
    else:
        ax.set_xlabel( r'$\sqrt{\lambda}$, $\omega$' )

    ax.set_ylabel( plot_labels['y_axis'] )

    if new_axes:
        ax.set_xlim( [fmin, fmax] )
        ax.set_ylim( plot_range )

    if show_legend:
        ax.legend()

    if show:
        pylab.show()
    return fig
    
def plot_gaps( fig_num, plot_rsc, gaps, kinds, freq_range,
               plot_range, show = False, clear = False, new_axes = False ):
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

def setup_band_gaps( pb, eigs, eig_vectors, opts, funmod ):
    """Setup density, volume info and eigenmomenta. Adjust frequency ranges."""
    dv_info = compute_density_volume_info( pb, opts.volume,
                                           opts.region_to_material )
    output( 'average density:', dv_info.average_density )
    
    if opts.fixed_eig_range is not None:
        mine, maxe = opts.fixed_eig_range
        ii = nm.where( (eigs > (mine**2.)) & (eigs < (maxe**2.)) )[0]
        freq_range_initial = nm.sqrt( eigs[ii] )
        opts.eig_range = (ii[0], ii[-1]+1) # +1 as it is a slice.
    else:
        freq_range_initial = nm.sqrt( eigs[slice( *opts.eig_range )] )
    output( 'initial freq. range     : [%8.3f, %8.3f]'\
            % tuple( freq_range_initial[[0,-1]] ) )

    if opts.eig_vector_transform is not None:
        fun = getattr( funmod, opts.eig_vector_transform[0] )
        def wrap_transform( vec, shape ):
            return fun( vec, shape, *opts.eig_vector_transform[1:] )
    else:
        wrap_transform = None
    output( 'mass matrix eigenmomenta...')
    pbar = MyBar( 'computing:' )
    tt = time.clock()
    aux = prepare_eigenmomenta( pb, opts.eigenmomentum,
                                opts.region_to_material,
                                eig_vectors, opts.teps, opts.teps_rel,
                                transform = wrap_transform, pbar = pbar )
    n_zeroed, valid, eigenmomenta = aux
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

    freq_info = Struct( name = 'freq_info',
                        freq_range_initial = freq_range_initial,
                        freq_range = freq_range,
                        freq_range_margins = freq_range_margins )

    return dv_info, eigenmomenta, n_zeroed, valid, freq_info
    
def detect_band_gaps( pb, eigs, eig_vectors, opts, funmod,
                      christoffel = None ):
    """Detect band gaps given solution to eigenproblem (eigs,
    eig_vectors). Only valid resonance frequencies (e.i. those for which
    corresponding eigenmomenta are above a given threshold) are taken into
    account.

    ??? make feps relative to ]f0, f1[ size ???
    """
    output( 'eigensolver:', opts.eigensolver )

    aux = setup_band_gaps( pb, eigs, eig_vectors, opts, funmod )
    dv_info, eigenmomenta, n_zeroed, valid, freq_info = aux

    fm = freq_info.freq_range_margins
    min_freq, max_freq = fm[0], fm[-1]
    output( 'freq. range with margins: [%8.3f, %8.3f]'\
            % (min_freq, max_freq) )

    df = opts.freq_step * (max_freq - min_freq)
    valid_eigenmomenta = eigenmomenta[valid,:]
    valid_eigs = eigs[valid]

    mass = AcousticMassTensor( valid_eigenmomenta, valid_eigs, dv_info )
    fz_callback = get_callback( mass, opts.eigensolver,
                                christoffel = christoffel, mode = 'find_zero' )
    trace_callback = get_callback( mass, opts.eigensolver,
                                   christoffel = christoffel, mode = 'trace' )

    n_col = 1 + (christoffel is not None)
    logs = [[] for ii in range( n_col + 1 )]
    gaps = []

    for ii in xrange( freq_info.freq_range.shape[0] + 1 ):

        f0, f1 = fm[[ii, ii+1]]
        output( 'interval: ]%.8f, %.8f[...' % (f0, f1) )

        f_delta = f1 - f0
        f_mid = 0.5 * (f0 + f1)
        if (f1 - f0) > (2.0 * opts.feps):
            num = min( 1000, max( 100, (f1 - f0) / df ) )
            a = nm.linspace( 0., 1., num )
            log_freqs = f0 + opts.feps \
                        + 0.5 * (nm.sin( (a - 0.5) * nm.pi ) + 1.0) \
                        * (f1 - f0 - 2.0 * opts.feps)
#            log_freqs = nm.linspace( f0 + opts.feps, f1 - opts.feps, num )
        else:
            log_freqs = nm.array( [f_mid - 1e-8 * f_delta,
                                   f_mid + 1e-8 * f_delta] )

        output( 'n_logged: %d' % log_freqs.shape[0] )

        log_mevp = [[] for ii in range( n_col )]
        for f in log_freqs:
            for ii, data in enumerate( trace_callback( f ) ):
                log_mevp[ii].append( data )

        # Get log for the first and last f in log_freqs.
        lf0 = log_freqs[0]
        lf1 = log_freqs[-1]

        log0, log1 = log_mevp[0][0], log_mevp[0][-1]
        min_eig0 = log0[0]
        max_eig1 = log1[-1]
        if min_eig0 > 0.0: # No gaps.
            gap = ([2, lf0, log0[0]], [2, lf0, log0[-1]])
        elif max_eig1 < 0.0: # Full interval strong gap.
            gap = ([1, lf1, log1[0]], [1, lf1, log1[-1]])
        else:
            llog_freqs = list( log_freqs )
            
            # Insert fmin, fmax into log.
            output( 'finding zero of the largest eig...' )
            smax, fmax, vmax = find_zero( lf0, lf1, fz_callback,
                                          opts.feps, opts.zeps, 1 )
            im = nm.searchsorted( log_freqs, fmax )
            llog_freqs.insert( im, fmax )
            for ii, data in enumerate( trace_callback( fmax ) ):
                log_mevp[ii].insert( im, data )

            output( '...done' )
            if smax in [0, 2]:
                output( 'finding zero of the smallest eig...' )
                # having fmax instead of f0 does not work if feps is large.
                smin, fmin, vmin = find_zero( lf0, lf1, fz_callback,
                                              opts.feps, opts.zeps, 0 )
                im = nm.searchsorted( log_freqs, fmin )
                # +1 due to fmax already inserted before.
                llog_freqs.insert( im+1, fmin )
                for ii, data in enumerate( trace_callback( fmin ) ):
                    log_mevp[ii].insert( im+1, data )

                output( '...done' )
            elif smax == 1:
                smin = 1 # both are negative everywhere.
                fmin, vmin = fmax, vmax

            gap = ([smin, fmin, vmin], [smax, fmax, vmax])

            log_freqs = nm.array( llog_freqs )

        output( gap[0] )
        output( gap[1] )
#        pause()
        gaps.append( gap )

        logs[0].append( log_freqs )
        for ii, data in enumerate( log_mevp ):
            logs[ii+1].append( nm.array( data, dtype = nm.float64 ) )

        output( '...done' )

    kinds = describe_gaps( gaps )

    return Struct( logs = logs, gaps = gaps, kinds = kinds,
                   valid = valid, eig_range = slice( *opts.eig_range ),
                   n_eigs = eigs.shape[0], n_zeroed = n_zeroed,
                   freq_range_initial = freq_info.freq_range_initial,
                   freq_range = freq_info.freq_range,
                   freq_range_margins = freq_info.freq_range_margins,
                   opts = opts )

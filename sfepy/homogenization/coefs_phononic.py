
import numpy as nm
import numpy.linalg as nla
import scipy as sc

from sfepy.base.base import output, get_default, dict_to_struct, assert_, Struct
from sfepy.base.timing import Timer
from sfepy.solvers import eig, Solver
from sfepy.linalg import norm_l2_along_axis
from sfepy.discrete.evaluate import eval_equations
from sfepy.homogenization.coefs_base import MiniAppBase, CorrMiniApp
from sfepy.homogenization.utils import coor_to_sym

def compute_eigenmomenta(em_equation, var_name, problem, eig_vectors,
                         transform=None):
    """
    Compute the eigenmomenta corresponding to given eigenvectors.
    """
    n_dof, n_eigs = eig_vectors.shape

    equations, variables = problem.create_evaluable(em_equation)
    variables.init_state()
    var = variables[var_name]

    n_c = var.n_components
    eigenmomenta = nm.empty((n_eigs, n_c), dtype=nm.float64)

    for ii in range(n_eigs):

        if transform is None:
            vec_phi, is_zero = eig_vectors[:,ii], False

        else:
            vec_phi, is_zero = transform(eig_vectors[:,ii], (n_dof / n_c, n_c))

        if is_zero:
            eigenmomenta[ii, :] = 0.0

        else:
            variables.set_state_parts({var_name: vec_phi})

            val = eval_equations(equations, variables)

            eigenmomenta[ii, :] = val

    return eigenmomenta

def get_ranges(freq_range, eigs):
    """
    Get an eigenvalue range slice and a corresponding initial frequency range
    within a given frequency range.
    """
    mine, maxe = freq_range
    ii = nm.where((eigs > (mine**2.)) & (eigs < (maxe**2.)))[0]
    freq_range_initial = nm.sqrt(eigs[ii])
    eig_range = (ii[0], ii[-1] + 1) # +1 as it is a slice.
    eig_range = slice(*eig_range)

    return freq_range_initial, eig_range

def cut_freq_range(freq_range, eigs, valid, freq_margins, eig_range,
                   fixed_freq_range, freq_eps):
    """
    Cut off masked resonance frequencies. Margins are preserved, like no
    resonances were cut.

    Returns
    -------
    freq_range : array
        The new range of frequencies.
    freq_range_margins : array
        The range of frequencies with prepended/appended margins equal to
        `fixed_freq_range` if it is not None.
    """
    n_eigs = eigs.shape[0]

    output('masked resonance frequencies in range:')
    output(nm.where(valid[eig_range] == False)[0])

    if fixed_freq_range is None:
        min_freq, max_freq = freq_range[0], freq_range[-1]
        margins = freq_margins * (max_freq - min_freq)
        prev_freq = min_freq - margins[0]
        next_freq = max_freq + margins[1]

        if eig_range.start > 0:
            prev_freq = max(nm.sqrt(eigs[eig_range.start - 1]) + freq_eps,
                            prev_freq)

        if eig_range.stop < n_eigs:
            next_freq = min(nm.sqrt(eigs[eig_range.stop]) - freq_eps,
                            next_freq)

        prev_freq = max(freq_eps, prev_freq)
        next_freq = max(freq_eps, next_freq, prev_freq + freq_eps)

    else:
        prev_freq, next_freq = fixed_freq_range

    freq_range = freq_range[valid[eig_range]]
    freq_range_margins = nm.r_[prev_freq, freq_range, next_freq]

    return freq_range, freq_range_margins

def split_chunks(indx):
    """Split index vector to chunks of consecutive numbers."""
    if not len(indx): return []

    delta = nm.ediff1d(indx, to_end=2)
    ir = nm.where(delta > 1)[0]

    chunks = []
    ic0 = 0
    for ic in ir:
        chunk = indx[ic0:ic+1]
        ic0 = ic + 1
        chunks.append(chunk)
    return chunks

def get_log_freqs(f0, f1, df, freq_eps, n_point_min, n_point_max):
    """
    Get logging frequencies.

    The frequencies get denser towards the interval boundaries.
    """
    f_delta = f1 - f0
    f_mid = 0.5 * (f0 + f1)

    if (f1 - f0) > (2.0 * freq_eps):
        num = min(n_point_max, max(n_point_min, int((f1 - f0) / df)))
        a = nm.linspace(0., 1., num)
        log_freqs = f0 + freq_eps \
                    + 0.5 * (nm.sin((a - 0.5) * nm.pi) + 1.0) \
                    * (f1 - f0 - 2.0 * freq_eps)

    else:
        log_freqs = nm.array([f_mid - 1e-8 * f_delta,
                              f_mid + 1e-8 * f_delta])

    return log_freqs

def detect_band_gaps(mass, freq_info, opts, gap_kind='normal', mtx_b=None):
    """
    Detect band gaps given solution to eigenproblem (eigs,
    eig_vectors). Only valid resonance frequencies (e.i. those for which
    corresponding eigenmomenta are above a given threshold) are taken into
    account.

    Notes
    -----
    - make freq_eps relative to ]f0, f1[ size?
    """
    output('eigensolver:', opts.eigensolver)

    fm = freq_info.freq_range_margins
    min_freq, max_freq = fm[0], fm[-1]
    output('freq. range with margins: [%8.3f, %8.3f]'
           % (min_freq, max_freq))

    df = opts.freq_step * (max_freq - min_freq)

    fz_callback = get_callback(mass.evaluate, opts.eigensolver,
                               mtx_b=mtx_b, mode='find_zero')
    trace_callback = get_callback(mass.evaluate, opts.eigensolver,
                                  mtx_b=mtx_b, mode='trace')

    n_col = 1 + (mtx_b is not None)
    logs = [[] for ii in range(n_col + 1)]
    gaps = []

    for ii in range(freq_info.freq_range.shape[0] + 1):

        f0, f1 = fm[[ii, ii+1]]
        output('interval: ]%.8f, %.8f[...' % (f0, f1))

        log_freqs = get_log_freqs(f0, f1, df, opts.freq_eps, 100, 1000)

        output('n_logged: %d' % log_freqs.shape[0])

        log_mevp = [[] for ii in range(n_col)]
        for f in log_freqs:
            for ii, data in enumerate(trace_callback(f)):
                log_mevp[ii].append(data)

        # Get log for the first and last f in log_freqs.
        lf0 = log_freqs[0]
        lf1 = log_freqs[-1]

        log0, log1 = log_mevp[0][0], log_mevp[0][-1]
        min_eig0 = log0[0]
        max_eig1 = log1[-1]
        if gap_kind == 'liquid':
            mevp = nm.array(log_mevp, dtype=nm.float64).squeeze()
            si = nm.where(mevp[:,0] < 0.0)[0]
            li = nm.where(mevp[:,-1] < 0.0)[0]
            wi = nm.setdiff1d(si, li)

            if si.shape[0] == 0: # No gaps.
                gap = ([2, lf0, log0[0]], [2, lf0, log0[-1]])
                gaps.append(gap)

            elif li.shape[0] == mevp.shape[0]: # Full interval strong gap.
                gap = ([1, lf1, log1[0]], [1, lf1, log1[-1]])
                gaps.append(gap)

            else:
                subgaps = []
                for chunk in split_chunks(li): # Strong gaps.
                    i0, i1 = chunk[0], chunk[-1]
                    fmin, fmax = log_freqs[i0], log_freqs[i1]
                    gap = ([1, fmin, mevp[i0,-1]], [1, fmax, mevp[i1,-1]])
                    subgaps.append(gap)

                for chunk in split_chunks(wi): # Weak gaps.
                    i0, i1 = chunk[0], chunk[-1]
                    fmin, fmax = log_freqs[i0], log_freqs[i1]
                    gap = ([0, fmin, mevp[i0,-1]], [2, fmax, mevp[i1,-1]])
                    subgaps.append(gap)
                gaps.append(subgaps)

        else:
            if min_eig0 > 0.0: # No gaps.
                gap = ([2, lf0, log0[0]], [2, lf0, log0[-1]])

            elif max_eig1 < 0.0: # Full interval strong gap.
                gap = ([1, lf1, log1[0]], [1, lf1, log1[-1]])

            else:
                llog_freqs = list(log_freqs)

                # Insert fmin, fmax into log.
                output('finding zero of the largest eig...')
                smax, fmax, vmax = find_zero(lf0, lf1, fz_callback,
                                             opts.freq_eps, opts.zero_eps, 1)
                im = nm.searchsorted(log_freqs, fmax)
                llog_freqs.insert(im, fmax)
                for ii, data in enumerate(trace_callback(fmax)):
                    log_mevp[ii].insert(im, data)

                output('...done')
                if smax in [0, 2]:
                    output('finding zero of the smallest eig...')
                    # having fmax instead of f0 does not work if freq_eps is
                    # large.
                    smin, fmin, vmin = find_zero(lf0, lf1, fz_callback,
                                                 opts.freq_eps, opts.zero_eps, 0)
                    im = nm.searchsorted(log_freqs, fmin)
                    # +1 due to fmax already inserted before.
                    llog_freqs.insert(im+1, fmin)
                    for ii, data in enumerate(trace_callback(fmin)):
                        log_mevp[ii].insert(im+1, data)

                    output('...done')

                elif smax == 1:
                    smin = 1 # both are negative everywhere.
                    fmin, vmin = fmax, vmax

                gap = ([smin, fmin, vmin], [smax, fmax, vmax])

                log_freqs = nm.array(llog_freqs)

            output(gap[0])
            output(gap[1])

            gaps.append(gap)

        logs[0].append(log_freqs)
        for ii, data in enumerate(log_mevp):
            logs[ii+1].append(nm.array(data, dtype = nm.float64))

        output('...done')

    kinds = describe_gaps(gaps)

    slogs = Struct(freqs=logs[0], eigs=logs[1])
    if n_col == 2:
        slogs.eig_vectors = logs[2]

    return slogs, gaps, kinds

def get_callback(mass, solver_kind, mtx_b=None, mode='trace'):
    r"""
    Return callback to solve band gaps or dispersion eigenproblem P.

    Notes
    -----
    Find zero callbacks return:
      eigenvalues

    Trace callbacks return:
      (eigenvalues,)
    or
      (eigenvalues, eigenvectors) (in full (dispoersion) mode)

    If `mtx_b` is None, the problem P is
      M w = \lambda w,
    otherwise it is
      omega^2 M w = \eta B w"""

    def find_zero_callback(f):
        meigs = eig(mass(f), eigenvectors=False, solver_kind=solver_kind)
        return meigs

    def find_zero_full_callback(f):
        meigs = eig((f**2) * mass(f), mtx_b=mtx_b,
                    eigenvectors=False, solver_kind=solver_kind)
        return meigs

    def trace_callback(f):
        meigs = eig(mass(f), eigenvectors=False, solver_kind=solver_kind)
        return meigs,

    def trace_full_callback(f):
        meigs, mvecs = eig((f**2) * mass(f), mtx_b=mtx_b,
                           eigenvectors=True, solver_kind=solver_kind)

        return meigs, mvecs

    if mtx_b is not None:
        mode += '_full'

    return eval(mode + '_callback')

def find_zero(f0, f1, callback, freq_eps, zero_eps, mode):
    r"""
    For f \in ]f0, f1[ find frequency f for which either the smallest (`mode` =
    0) or the largest (`mode` = 1) eigenvalue of problem P given by `callback`
    is zero.

    Returns
    -------
    flag : 0, 1, or 2
        The flag, see Notes below.
    frequency : float
        The found frequency.
    eigenvalue : float
        The eigenvalue corresponding to the found frequency.

    Notes
    -----
    Meaning of the return value combinations:

    =====  ======  ========
    mode    flag    meaning
    =====  ======  ========
    0, 1    0       eigenvalue -> 0 for f \in ]f0, f1[
    0       1       f -> f1, smallest eigenvalue < 0
    0       2       f -> f0, smallest eigenvalue > 0 and -> -\infty
    1       1       f -> f1, largest eigenvalue < 0 and  -> +\infty
    1       2       f -> f0, largest eigenvalue > 0
    =====  ======  ========
    """
    fm, fp = f0, f1
    ieig = {0 : 0, 1 : -1}[mode]
    while 1:
        f = 0.5 * (fm + fp)
        meigs = callback(f)

        val = meigs[ieig]
        ## print f, f0, f1, fm, fp, val
        ## print '%.16e' % f, '%.16e' % fm, '%.16e' % fp, '%.16e' % val

        if ((abs(val) < zero_eps)
            or ((fp - fm) < (abs(fm) * nm.finfo(float).eps))):
            return 0, f, val

        if mode == 0:
            if (f - f0) < freq_eps:
                return 2, f0, val

            elif (f1 - f) < freq_eps:
                return 1, f1, val

        elif mode == 1:
            if (f1 - f) < freq_eps:
                return 1, f1, val

            elif (f - f0) < freq_eps:
                return 2, f0, val

        if val > 0.0:
            fp = f

        else:
            fm = f

def describe_gaps(gaps):
    kinds = []
    for ii, gap in enumerate(gaps):
        if isinstance(gap, list):
            subkinds = []
            for gmin, gmax in gap:
                if (gmin[0] == 2) and (gmax[0] == 2):
                    kind = ('p', 'propagation zone')
                elif (gmin[0] == 1) and (gmax[0] == 1):
                    kind = ('is', 'inner strong band gap')
                elif (gmin[0] == 0) and (gmax[0] == 2):
                    kind = ('iw', 'inner weak band gap')
                subkinds.append(kind)
            kinds.append(subkinds)

        else:
            gmin, gmax = gap

            if (gmin[0] == 2) and (gmax[0] == 2):
                kind = ('p', 'propagation zone')
            elif (gmin[0] == 1) and (gmax[0] == 2):
                kind = ('w', 'full weak band gap')
            elif (gmin[0] == 0) and (gmax[0] == 2):
                kind = ('wp', 'weak band gap + propagation zone')
            elif (gmin[0] == 1) and (gmax[0] == 1):
                kind = ('s', 'full strong band gap (due to end of freq.'
                        ' range or too large thresholds)')
            elif (gmin[0] == 1) and (gmax[0] == 0):
                kind = ('sw', 'strong band gap + weak band gap')
            elif (gmin[0] == 0) and (gmax[0] == 0):
                kind = ('swp', 'strong band gap + weak band gap +'
                        ' propagation zone')
            else:
                msg = 'impossible band gap combination: %d, %d' % (gmin, gmax)
                raise ValueError(msg)
            kinds.append(kind)

    return kinds

def get_gap_ranges(freq_range, gaps, kinds):
    """
    For each (potential) band gap in `gaps`, return the frequency ranges of its
    parts according to `kinds`.
    """

    def get_ranges(ii, f0, f1, kind, kind_desc, gmin, gmax):
        if kind == 'p':
            ranges = [(f0, f1)]
        elif kind == 'w':
            ranges = [(f0, f1)]
        elif kind == 'wp':
            ranges = [(f0, gmin[1]), (gmin[1], f1)]
        elif kind == 's':
            ranges = [(f0, f1)]
        elif kind == 'sw':
            ranges = [(f0, gmax[1]), (gmax[1], f1)]
        elif kind == 'swp':
            ranges = [(f0, gmax[1]), (gmax[1], gmin[1]), (gmin[1], f1)]
        elif kind == 'is':
            ranges = [(gmin[1], gmax[1])]
        elif kind == 'iw':
            ranges = [(gmin[1], gmax[1])]
        else:
            msg = 'impossible band gap combination! (%f %d)' % (gmin, gmax)
            raise ValueError(msg)

        return ranges

    gap_ranges = []
    for ii in range(len(freq_range) - 1):
        f0, f1 = freq_range[[ii, ii+1]]
        gap = gaps[ii]

        if isinstance(gap, list):
            ranges = []
            for ig, (gmin, gmax) in enumerate(gap):
                kind, kind_desc = kinds[ii][ig]
                aux = get_ranges(ii, f0, f1, kind, kind_desc, gmin, gmax)
                ranges.append(aux)

        else:
            gmin, gmax = gap
            kind, kind_desc = kinds[ii]
            ranges = get_ranges(ii, f0, f1, kind, kind_desc, gmin, gmax)

        gap_ranges.append(ranges)

    return gap_ranges

def compute_cat_sym_sym(coef, iw_dir):
    """
    Christoffel acoustic tensor (part) of elasticity tensor dimension.
    """
    dim = iw_dir.shape[0]

    cat = nm.zeros((dim, dim), dtype=nm.float64)
    for ii in range(dim):
        for ij in range(dim):
            ir = coor_to_sym(ii, ij, dim)
            for ik in range(dim):
                for il in range(dim):
                    ic = coor_to_sym(ik, il, dim)
                    cat[ii,ik] += coef[ir,ic] * iw_dir[ij] * iw_dir[il]

    return cat

def compute_cat_dim_sym(coef, iw_dir):
    """
    Christoffel acoustic tensor part of piezo-coupling tensor dimension.
    """
    dim = iw_dir.shape[0]

    cat = nm.zeros((dim,), dtype=nm.float64)
    for ii in range(dim):
        for ij in range(dim):
            ir = coor_to_sym(ii, ij, dim)
            for ik in range(dim):
                cat[ii] += coef[ik,ir] * iw_dir[ij] * iw_dir[ik]

    return cat

def compute_cat_dim_dim(coef, iw_dir):
    """
    Christoffel acoustic tensor part of dielectric tensor dimension.
    """
    cat = nm.dot(nm.dot(coef, iw_dir), iw_dir)

    return cat

class SimpleEVP(CorrMiniApp):
    """
    Simple eigenvalue problem.
    """

    def process_options(self):
        get = self.options.get

        return Struct(eigensolver=get('eigensolver', 'eig.sgscipy'),
                      elasticity_contrast=get('elasticity_contrast', 1.0),
                      scale_epsilon=get('scale_epsilon', 1.0),
                      save_eig_vectors=get('save_eig_vectors', (0, 0)))

    def __call__(self, problem=None, data=None):
        problem = get_default(problem, self.problem)
        opts = self.app_options

        if self.equations is not None:
            problem.set_equations(self.equations)
            problem.select_bcs(ebc_names=self.ebcs, epbc_names=self.epbcs,
                            lcbc_names=self.get('lcbcs', []))
            problem.update_materials(problem.ts)

            self.init_solvers(problem)

        mtx_a, mtx_m, data = self.prepare_matrices(problem)

        output('computing resonance frequencies...')
        tt = [0]

        if isinstance(mtx_a, sc.sparse.spmatrix):
            mtx_a = mtx_a.toarray()
        if isinstance(mtx_m, sc.sparse.spmatrix):
            mtx_m = mtx_m.toarray()

        eigs, mtx_s_phi = eig(mtx_a, mtx_m, return_time=tt,
                              solver_kind=opts.eigensolver)
        eigs[eigs<0.0] = 0.0
        output('...done in %.2f s' % tt[0])
        output('original eigenfrequencies:')
        output(eigs)
        opts = self.app_options
        epsilon2 = opts.scale_epsilon * opts.scale_epsilon
        eigs_rescaled = (opts.elasticity_contrast / epsilon2)  * eigs
        output('rescaled eigenfrequencies:')
        output(eigs_rescaled)
        output('number of eigenfrequencies: %d' % eigs.shape[0])

        try:
            assert_(nm.isfinite(eigs).all())
        except ValueError:
            from sfepy.base.base import debug; debug()

        mtx_phi, eig_vectors = self.post_process(eigs, mtx_s_phi, data,
                                                 problem)

        self.save(eigs, mtx_phi, problem)

        evp = Struct(name='evp', eigs=eigs, eigs_rescaled=eigs_rescaled,
                     eig_vectors=eig_vectors)

        return evp

    def prepare_matrices(self, problem):
        mtx_a = problem.evaluate(self.equations['lhs'], mode='weak',
                                 auto_init=True, dw_mode='matrix')

        mtx_m = problem.evaluate(self.equations['rhs'], mode='weak',
                                 dw_mode='matrix')

        return mtx_a, mtx_m, None

    def post_process(self, eigs, mtx_s_phi, data, problem):
        n_eigs = eigs.shape[0]

        variables = problem.get_variables()

        mtx_phi = nm.empty((variables.di.n_dof_total, mtx_s_phi.shape[1]),
                           dtype=nm.float64)

        make_full = variables.make_full_vec
        for ii in range(n_eigs):
            mtx_phi[:,ii] = make_full(mtx_s_phi[:,ii])

        return mtx_phi, mtx_phi

    def save(self, eigs, mtx_phi, problem):
        save = self.app_options.save_eig_vectors

        n_eigs = eigs.shape[0]

        out = {}
        variables = problem.set_default_state()
        for ii in range(n_eigs):
            if (ii >= save[0]) and (ii < (n_eigs - save[1])): continue
            variables.set_state(mtx_phi[:,ii], force=True)
            aux = variables.create_output()
            for name, val in aux.items():
                out[name+'%03d' % ii] = val

        if self.post_process_hook is not None:
            out = self.post_process_hook(out, problem, mtx_phi)

        problem.domain.mesh.write(self.save_name + '.vtk', io='auto', out=out)

        fd = open(self.save_name + '_eigs.txt', 'w')
        eigs.tofile(fd, ' ')
        fd.close()

class SchurEVP(SimpleEVP):
    """
    Schur complement eigenvalue problem.
    """

    def prepare_matrices(self, problem):
        """
        A = K + B^T D^{-1} B
        """
        equations = problem.equations
        mtx = equations.eval_tangent_matrices(None, problem.mtx_a,
                                              by_blocks=True)

        ls = Solver.any_from_conf(problem.ls_conf
                                  + Struct(use_presolve=True), mtx=mtx['D'])

        mtx_b, mtx_m = mtx['B'], mtx['M']
        mtx_dib = nm.empty(mtx_b.shape, dtype=mtx_b.dtype)
        for ic in range(mtx_b.shape[1]):
            mtx_dib[:,ic] = ls(mtx_b[:,ic].toarray().squeeze())
        mtx_a = mtx['K'] + mtx_b.T * mtx_dib

        return mtx_a, mtx_m, mtx_dib

    def post_process(self, eigs, mtx_s_phi, mtx_dib, problem):
        n_eigs = eigs.shape[0]

        variables = problem.get_variables()

        mtx_phi = nm.empty((variables.di.n_dof_total, mtx_s_phi.shape[1]),
                           dtype=nm.float64)

        make_full = variables.make_full_vec

        # Update also eliminated variables.
        schur = self.app_options.schur
        primary_var = schur['primary_var']
        eliminated_var = schur['eliminated_var']

        mtx_s_phi_schur = - sc.dot(mtx_dib, mtx_s_phi)
        aux = nm.empty((variables.adi.n_dof_total,), dtype=nm.float64)
        setv = variables.set_vec_part
        for ii in range(n_eigs):
            setv(aux, primary_var, mtx_s_phi[:,ii], reduced=True)
            setv(aux, eliminated_var, mtx_s_phi_schur[:,ii],
                 reduced=True)

            mtx_phi[:,ii] = make_full(aux)

        indx = variables.get_indx(primary_var)
        eig_vectors = mtx_phi[indx,:]

        return mtx_phi, eig_vectors

class DensityVolumeInfo(MiniAppBase):
    """
    Determine densities of regions specified in `region_to_material`, and
    compute average density based on region volumes.
    """

    def __call__(self, volume=None, problem=None, data=None):
        problem = get_default(problem, self.problem)

        vf = data[self.requires[0]]

        average_density = 0.0
        total_volume = 0.0
        volumes = {}
        densities = {}
        for region_name, aux in self.region_to_material.items():
            vol = vf['volume_' + region_name]

            mat_name, item_name = aux
            conf = problem.conf.get_item_by_name('materials', mat_name)
            density = conf.values[item_name]

            output('region %s: volume %f, density %f' % (region_name,
                                                         vol, density))

            volumes[region_name] = vol
            densities[region_name] = density

            average_density += vol * density
            total_volume += vol

        true_volume = self._get_volume(volume)

        output('total volume:', true_volume)

        average_density /= true_volume

        return Struct(name='density_volume_info',
                      average_density=average_density,
                      total_volume=true_volume,
                      volumes=volumes,
                      densities=densities,
                      to_file_txt=self.to_file_txt)

    @staticmethod
    def to_file_txt(fd, float_format, dv_info):
        ff = float_format + '\n'

        fd.write('total volume:\n')
        fd.write(ff % dv_info.total_volume)
        fd.write('average density:\n')
        fd.write(ff % dv_info.average_density)

        for key, val in dv_info.volumes.items():
            fd.write('%s volume:\n' % key)
            fd.write(ff % val)
            fd.write('%s density:\n' % key)
            fd.write(ff % dv_info.densities[key])

class Eigenmomenta(MiniAppBase):
    """
    Eigenmomenta corresponding to eigenvectors.

    Parameters
    ----------
    var_name : str
        The name of the variable used in the integral.
    threshold : float
        The threshold under which an eigenmomentum is considered zero.
    threshold_is_relative : bool
        If True, the `threshold` is relative w.r.t. max. norm of eigenmomenta.
    transform : callable, optional
        Optional function for transforming the eigenvectors before computing
        the eigenmomenta.

    Returns
    -------
    eigenmomenta : Struct
        The resulting eigenmomenta. An eigenmomentum above threshold is marked
        by the attribute 'valid' set to True.
    """

    def process_options(self):
        options = dict_to_struct(self.options)
        get = options.get

        return Struct(var_name=get('var_name', None,
                                   'missing "var_name" in options!'),
                      threshold=get('threshold', 1e-4),
                      threshold_is_relative=get('threshold_is_relative', True),
                      transform=get('transform', None))

    def __call__(self, volume=None, problem=None, data=None):
        problem = get_default(problem, self.problem)
        opts = self.app_options

        evp, dv_info = [data[ii] for ii in self.requires]

        output('computing eigenmomenta...')

        if opts.transform is not None:
            fun = problem.conf.get_function(opts.transform[0])
            def wrap_transform(vec, shape):
                return fun(vec, shape, *opts.eig_vector_transform[1:])

        else:
            wrap_transform = None

        timer = Timer(start=True)
        eigenmomenta = compute_eigenmomenta(self.expression, opts.var_name,
                                            problem, evp.eig_vectors,
                                            wrap_transform)
        output('...done in %.2f s' % timer.stop())

        n_eigs = evp.eigs.shape[0]

        mag = norm_l2_along_axis(eigenmomenta)

        if opts.threshold_is_relative:
            tol = opts.threshold * mag.max()
        else:
            tol = opts.threshold

        valid = nm.where(mag < tol, False, True)
        mask = nm.where(valid == False)[0]
        eigenmomenta[mask, :] = 0.0
        n_zeroed = mask.shape[0]

        output('%d of %d eigenmomenta zeroed (under %.2e)'\
                % (n_zeroed, n_eigs, tol))

        out = Struct(name='eigenmomenta', n_zeroed=n_zeroed,
                     eigenmomenta=eigenmomenta, valid=valid,
                     to_file_txt=None)
        return out

class AcousticMassTensor(MiniAppBase):
    """
    The acoustic mass tensor for a given frequency.

    Returns
    -------
    self : AcousticMassTensor instance
        This class instance whose `evaluate()` method computes for a given
        frequency the required tensor.

    Notes
    -----
    `eigenmomenta`, `eigs` should contain only valid resonances.
    """
    to_file_txt = None

    def __call__(self, volume=None, problem=None, data=None):
        evp, self.dv_info, ema = [data[ii] for ii in self.requires]

        self.eigs = evp.eigs[ema.valid]
        self.eigenmomenta = ema.eigenmomenta[ema.valid, :]

        return self

    def evaluate(self, freq):
        ema = self.eigenmomenta

        n_c = ema.shape[1]
        fmass = nm.zeros((n_c, n_c), dtype=nm.float64)

        num, denom = self.get_coefs(freq)
        de = 1.0 / denom
        if not nm.isfinite(de).all():
            raise ValueError('frequency %e too close to resonance!' % freq)

        for ir in range(n_c):
            for ic in range(n_c):
                if ir <= ic:
                    val = nm.sum(num * de * (ema[:, ir] * ema[:, ic]))
                    fmass[ir, ic] += val
                else:
                    fmass[ir, ic] = fmass[ic, ir]

        eye = nm.eye(n_c, n_c, dtype=nm.float64)
        mtx_mass = (eye * self.dv_info.average_density) \
                   - (fmass / self.dv_info.total_volume)

        return mtx_mass

    def get_coefs(self, freq):
        """
        Get frequency-dependent coefficients.
        """
        f2 = freq*freq
        de = f2 - self.eigs
        return f2, de

class AcousticMassLiquidTensor(AcousticMassTensor):

    def get_coefs(self, freq):
        """
        Get frequency-dependent coefficients.
        """
        eigs = self.eigs

        f2 = freq*freq
        aux = (f2 - self.gamma * eigs)
        num = f2 * aux
        denom = aux*aux + f2*(self.eta*self.eta)*nm.power(eigs, 2.0)
        return num, denom

class AppliedLoadTensor(AcousticMassTensor):
    """
    The applied load tensor for a given frequency.

    Returns
    -------
    self : AppliedLoadTensor instance
        This class instance whose `evaluate()` method computes for a given
        frequency the required tensor.

    Notes
    -----
    `eigenmomenta`, `ueigenmomenta`, `eigs` should contain only valid
    resonances.
    """
    to_file_txt = None

    def __call__(self, volume=None, problem=None, data=None):
        evp, self.dv_info, ema, uema = [data[ii] for ii in self.requires]

        self.eigs = evp.eigs[ema.valid]
        self.eigenmomenta = ema.eigenmomenta[ema.valid, :]
        self.ueigenmomenta = uema.eigenmomenta[uema.valid, :]

        return self

    def evaluate(self, freq):
        ema, uema = self.eigenmomenta, self.ueigenmomenta

        n_c = ema.shape[1]
        fload = nm.zeros((n_c, n_c), dtype=nm.float64)

        num, denom = self.get_coefs(freq)
        de = 1.0 / denom
        if not nm.isfinite(de).all():
            raise ValueError('frequency %e too close to resonance!' % freq)

        for ir in range(n_c):
            for ic in range(n_c):
                val = nm.sum(num * de * (ema[:, ir] * uema[:, ic]))
                fload[ir, ic] += val

        eye = nm.eye(n_c, n_c, dtype=nm.float64)

        mtx_load = eye - (fload / self.dv_info.total_volume)

        return mtx_load

class BandGaps(MiniAppBase):
    """
    Band gaps detection.

    Parameters
    ----------
    eigensolver : str
        The name of the eigensolver for mass matrix eigenvalues.
    eig_range : (int, int)
        The eigenvalues range (squared frequency) to consider.
    freq_margins : (float, float)
        Margins in percents of initial frequency range given by
        `eig_range` by which the range is increased.
    fixed_freq_range : (float, float)
        The frequency range to consider. Has precedence over `eig_range`
        and `freq_margins`.
    freq_step : float
        The frequency step for tracing, in percent of the frequency range.
    freq_eps : float
        The frequency difference smaller than `freq_eps` is considered zero.
    zero_eps : float
        The tolerance for finding zeros of mass matrix eigenvalues.
    detect_fun : callable
        The function for detecting the band gaps. Default is
        :func:`detect_band_gaps()`.
    log_save_name : str
        If not None, the band gaps log is to be saved under the given name.
    raw_log_save_name : str
        If not None, the raw band gaps log is to be saved under the given name.
    """

    def process_options(self):
        get = self.options.get

        freq_margins = get('freq_margins', (5, 5))
        # Given per cent.
        freq_margins = 0.01 * nm.array(freq_margins, dtype=nm.float64)

        # Given in per cent.
        freq_step = 0.01 * get('freq_step', 5)

        return Struct(eigensolver=get('eigensolver', 'eig.sgscipy'),
                      eig_range=get('eig_range', None),
                      freq_margins=freq_margins,
                      fixed_freq_range=get('fixed_freq_range', None),
                      freq_step=freq_step,

                      freq_eps=get('freq_eps', 1e-8),
                      zero_eps=get('zero_eps', 1e-8),
                      detect_fun=get('detect_fun', detect_band_gaps),
                      log_save_name=get('log_save_name', None),
                      raw_log_save_name=get('raw_log_save_name', None))

    def __call__(self, volume=None, problem=None, data=None):
        problem = get_default(problem, self.problem)
        opts = self.app_options

        evp, ema, mass = [data[ii] for ii in self.requires[:3]]
        if len(self.requires) == 4:
            mtx_b = data[self.requires[3]]

        else:
            mtx_b = None

        eigs = evp.eigs

        self.fix_eig_range(eigs.shape[0])

        if opts.fixed_freq_range is not None:
            (freq_range_initial,
             opts.eig_range) = get_ranges(opts.fixed_freq_range, eigs)

        else:
            opts.eig_range = slice(*opts.eig_range)
            freq_range_initial = nm.sqrt(eigs[opts.eig_range])

        output('initial freq. range     : [%8.3f, %8.3f]'
               % tuple(freq_range_initial[[0, -1]]))

        aux = cut_freq_range(freq_range_initial, eigs, ema.valid,
                             opts.freq_margins, opts.eig_range,
                             opts.fixed_freq_range,
                             opts.freq_eps)
        freq_range, freq_range_margins = aux
        if len(freq_range):
            output('freq. range             : [%8.3f, %8.3f]'
                   % tuple(freq_range[[0, -1]]))

        else:
            # All masked.
            output('freq. range             : all masked!')

        freq_info = Struct(name='freq_info',
                           freq_range_initial=freq_range_initial,
                           freq_range=freq_range,
                           freq_range_margins=freq_range_margins)

        logs, gaps, kinds = opts.detect_fun(mass, freq_info, opts, mtx_b=mtx_b)
        gap_ranges = get_gap_ranges(freq_range_margins, gaps, kinds)

        bg = Struct(name='band_gaps', logs=logs, gaps=gaps, kinds=kinds,
                    gap_ranges=gap_ranges,
                    valid=ema.valid, eig_range=opts.eig_range,
                    n_eigs=eigs.shape[0], n_zeroed=ema.n_zeroed,
                    freq_range_initial=freq_info.freq_range_initial,
                    freq_range=freq_info.freq_range,
                    freq_range_margins=freq_info.freq_range_margins,
                    opts=opts, to_file_txt=self.to_file_txt,
                    log_save_name=opts.log_save_name,
                    raw_log_save_name=opts.raw_log_save_name,
                    save_log=self.save_log)

        return bg

    def fix_eig_range(self, n_eigs):
        eig_range = get_default(self.app_options.eig_range, (0, n_eigs))
        if eig_range[-1] < 0:
            eig_range[-1] += n_eigs + 1

        assert_(eig_range[0] < (eig_range[1] - 1))
        assert_(eig_range[1] <= n_eigs)
        self.app_options.eig_range = eig_range

    @staticmethod
    def to_file_txt(fd, float_format, bg):
        if bg.log_save_name is not None:
            fd.write(bg.log_save_name + '\n')

        else:
            fd.write('--\n')

    @staticmethod
    def save_log(filename, float_format, bg):
        """
        Save band gaps, valid flags and eigenfrequencies.
        """
        fd = open(filename, 'w')
        freq_range = bg.freq_range_margins
        fd.write('n_zeroed: %d\n' % bg.n_zeroed)
        fd.write('n_eigs: %d\n' % bg.n_eigs)
        n_row = len(freq_range) - 1
        fd.write('n_ranges: %d\n' % n_row)
        fd.write('f0 f1 flag_min f_min v_min flag_max f_max v_max'
                  ' kind\ndesc\n')

        ff = float_format
        format = "%s %s %%d %s %s %%d %s %s %%s\n%%s\n" % (6 * (ff,))
        for ir in range(n_row):
            f0, f1 = freq_range[[ir, ir+1]]
            gmin, gmax = bg.gaps[ir]
            fd.write(format % ((f0, f1) + tuple(gmin) + tuple(gmax)
                                + bg.kinds[ir]))

        fd.write('\nname kind f_from f_to index f0 f1\n')
        format = '%%s %%s %s %s %%d %s %s\n' % (4 * (ff,))
        for ii, f0 in enumerate(bg.freq_range_margins[:-1]):
            f1 = bg.freq_range_margins[ii + 1]
            kind = bg.kinds[ii][0]
            for ir, rng in enumerate(bg.gap_ranges[ii]):
                fd.write(format
                         % (bg.name, kind[ir], rng[0], rng[1], ii, f0, f1))

        n_row = len(freq_range)
        fd.write('\nn_resonance: %d\n' % n_row)
        fd.write('valid f\n')
        freq_range = bg.freq_range_initial
        valid_in_range = bg.valid[bg.eig_range]
        format = "%%d %s\n" % ff
        for ir in range(n_row):
            fd.write(format % (valid_in_range[ir], freq_range[ir]))
        fd.close()

class ChristoffelAcousticTensor(MiniAppBase):

    def process_options(self):
        get = self.options.get
        return Struct(mode=get('mode', 'simple'),
                      incident_wave_dir=get('incident_wave_dir', None))

    r"""
    Compute Christoffel acoustic tensor (cat) given the incident wave
    direction (unit vector).

    Parameters
    ----------
    mode : 'simple' or 'piezo'
        The call mode.
    incident_wave_dir : array
        The incident wave direction vector.

    Returns
    -------
    cat : array
        The Christoffel acoustic tensor.

    Notes
    -----
    - If mode == 'simple', only the elasticity tensor :math:`C_{ijkl}` is used
      and cat := :math:`\Gamma_{ik} = C_{ijkl} n_j n_l`.

    - If mode == 'piezo', also the piezo-coupling :math:`G_{ijk}` and
      dielectric :math:`D_{ij}` tensors are used and cat := :math:`H_{ik} =
      \Gamma_{ik} + \frac{1}{\xi} \gamma_i \gamma_j`, where :math:`\gamma_i =
      G_{kij} n_j n_k`, :math:`\xi = D_{kl} n_k n_l`.
    """

    def __call__(self, volume=None, problem=None, data=None):
        problem = get_default(problem, self.problem)
        opts = self.app_options

        iw_dir = nm.array(opts.incident_wave_dir, dtype=nm.float64)
        dim = problem.get_dim()
        assert_(dim == iw_dir.shape[0])

        iw_dir = iw_dir / nla.norm(iw_dir)

        elastic = data[self.requires[0]]

        cat = compute_cat_sym_sym(elastic, iw_dir)

        if opts.mode =='piezo':
            dielectric, coupling = [data[ii] for ii in self.requires[1:]]
            xi = compute_cat_dim_dim(dielectric, iw_dir)
            gamma = compute_cat_dim_sym(coupling, iw_dir)
            cat += nm.outer(gamma, gamma) / xi

        return cat

class PolarizationAngles(MiniAppBase):
    """
    Compute polarization angles, i.e., angles between incident wave direction
    and wave vectors. Vector length does not matter - eigenvectors are used
    directly.
    """

    def process_options(self):
        get = self.options.get
        return Struct(incident_wave_dir=get('incident_wave_dir', None))

    def __call__(self, volume=None, problem=None, data=None):
        problem = get_default(problem, self.problem)
        opts = self.app_options

        iw_dir = nm.array(opts.incident_wave_dir, dtype=nm.float64)
        dim = problem.get_dim()
        assert_(dim == iw_dir.shape[0])

        iw_dir = iw_dir / nla.norm(iw_dir)

        dispersion = data[self.requires[0]]

        wave_vectors = dispersion.logs.eig_vectors

        pas = []

        iw_dir = iw_dir / nla.norm(iw_dir)
        idims = list(range(iw_dir.shape[0]))
        pi2 = 0.5 * nm.pi
        for vecs in wave_vectors:
            pa = nm.empty(vecs.shape[:-1], dtype=nm.float64)
            for ir, vec in enumerate(vecs):
                for ic in idims:
                    vv = vec[:,ic]
                    # Ensure the angle is in [0, pi/2].
                    val = nm.arccos(nm.dot(iw_dir, vv) / nla.norm(vv))
                    if val > pi2:
                        val = nm.pi - val
                    pa[ir,ic] = val

            pas.append(pa)

        return pas

class PhaseVelocity(MiniAppBase):
    """
    Compute phase velocity.
    """

    def process_options(self):
        get = self.options.get

        return Struct(eigensolver=get('eigensolver', 'eig.sgscipy'))

    def __call__(self, volume=None, problem=None, data=None):
        problem = get_default(problem, self.problem)
        opts = self.app_options

        dv_info, cat = [data[ii] for ii in self.requires]

        output('average density:', dv_info.average_density)

        dim = problem.get_dim()
        eye = nm.eye(dim, dim, dtype=nm.float64)
        mtx_mass = eye * dv_info.average_density

        meigs, mvecs = eig(mtx_mass, mtx_b=cat,
                           eigenvectors=True, solver_kind=opts.eigensolver)
        phase_velocity = 1.0 / nm.sqrt(meigs)

        return phase_velocity

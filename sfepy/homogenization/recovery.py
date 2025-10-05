import os

import numpy as nm

from sfepy.base.base import get_default, Struct, output
from sfepy.base.ioutils import get_print_info
from sfepy.base.timing import Timer
from sfepy.discrete.fem import extend_cell_data, Mesh
from sfepy.homogenization.utils import coor_to_sym
from sfepy.base.conf import get_standard_keywords
from sfepy.discrete import Problem, Region
from sfepy.base.conf import ProblemConf
from sfepy.homogenization.coefficients import Coefficients
from sfepy.homogenization.micmac import get_correctors_from_file_hdf5
import os.path as op
import atexit

shared = Struct()

#
# TODO : interpolate fvars to macro times. ?mid-points?
#
# TODO : clean-up!
#


def get_output_suffix(iel, ts, naming_scheme, format, output_format):
    if output_format != 'h5':
        if naming_scheme == 'step_iel':
            suffix = '.'.join((ts.suffix % ts.step, format % iel))
        else:
            suffix = '.'.join((format % iel, ts.suffix % ts.step))

    else:
        suffix = format % iel

    return suffix


def convolve_field_scalar(fvars, pvars, iel, ts):
    r"""
    .. math::
      \int_0^t f(t-s) p(s) ds

    Notes
    -----
    - t is given by step
    - f: fvars
      scalar field variables, defined in a micro domain, have shape [step][fmf
      dims]
    - p: pvars
      scalar point variables, a scalar in a point of macro-domain, FMField
      style have shape [n_step][var dims]
    """

    step0 = max(0, ts.step - fvars.steps[-1])

    val = nm.zeros_like(fvars[0])
    for ik in range(step0, ts.step + 1):
        vf = fvars[ts.step-ik]
        vp = pvars[ik][iel, 0, 0, 0]
        val += vf * vp * ts.dt

    return val


def convolve_field_sym_tensor(fvars, pvars, var_name, dim, iel, ts):
    r"""
    .. math::
      \int_0^t f^{ij}(t-s) p_{ij}(s) ds

    Notes
    -----
    - t is given by step
    - f: fvars
      field variables, defined in a micro domain, have shape [step][fmf dims]
    - p: pvars
      sym. tensor point variables, a scalar in a point of
      macro-domain, FMField style, have shape [dim, dim][var_name][n_step][var
      dims]
    """

    step0 = max(0, ts.step - fvars[0, 0][var_name].steps[-1])

    val = nm.zeros_like(fvars[0, 0][var_name][0])
    for ik in range(step0, ts.step + 1):
        for ir in range(dim):
            for ic in range(dim):
                ii = coor_to_sym(ir, ic, dim)
                vf = fvars[ir, ic][var_name][ts.step-ik]
                vp = pvars[ik][iel, 0, ii, 0]
                val += vf * vp * ts.dt
    return val


def add_strain_rs(corrs_rs, strain, vu, dim, iel, out=None):
    if out is None:
        out = nm.zeros_like(corrs_rs[0, 0][vu][0])

    for ir in range(dim):
        for ic in range(dim):
            ii = coor_to_sym(ir, ic, dim)
            out += corrs_rs[ir, ic][vu].data * strain[iel, 0, ii, 0]
    return out


def combine_scalar_grad(corrs, grad, vn, ii, shift_coors=None):
    r"""
    .. math::
      \eta_k \partial_k^x p

    or

    .. math::
      (y_k + \eta_k) \partial_k^x p
    """
    dim = grad.shape[2]

    if shift_coors is None:
        out = corrs[0][vn].data * grad[ii, 0, 0, 0]
        for ir in range(1, dim):
            out += corrs[ir][vn].data * grad[ii, 0, ir, 0]

    else:
        out = (shift_coors[:, 0] + corrs[0][vn].data) * grad[ii, 0, 0, 0]
        for ir in range(1, dim):
            out += (shift_coors[:, ir] + corrs[ir][vn].data) \
                * grad[ii, 0, ir, 0]

    return out


def compute_u_corr_steady(corrs_rs, strain, vu, dim, iel):
    r"""
    .. math::
      \sum_{ij} \bm{\omega}^{ij}\, e_{ij}(\bm{u})

    Notes
    -----
    - iel = element number
    """
    u_corr = add_strain_rs(corrs_rs, strain, vu, dim, iel)
    return u_corr


def compute_u_corr_time(corrs_rs, dstrains, corrs_pressure, pressures,
                        vu, dim, iel, ts):
    r"""
    .. math::
      \sum_{ij} \left[ \int_0^t \bm{\omega}^{ij}(t-s) {\mathrm{d} \over
      \mathrm{d} s} e_{ij}(\bm{u}(s))\,ds\right] + \int_0^t
      \widetilde{\bm{\omega}}^P(t-s)\,p(s)\,ds
    """
    u_corr = convolve_field_scalar(corrs_pressure[vu], pressures,
                                   iel, ts)
    u_corr += convolve_field_sym_tensor(corrs_rs, dstrains, vu,
                                        dim, iel, ts)
    return u_corr


def compute_p_corr_steady(corrs_pressure, pressure, vp, iel):
    r"""
    .. math::
      \widetilde\pi^P\,p
    """
    p_corr = corrs_pressure[vp].data * pressure[iel, 0, 0, 0]
    return p_corr


def compute_p_corr_time(corrs_rs, dstrains, corrs_pressure, pressures,
                        vdp, dim, iel, ts):
    r"""
    .. math::
      \sum_{ij} \int_0^t {\mathrm{d} \over \mathrm{d} t}
      \widetilde\pi^{ij}(t-s)\, {\mathrm{d} \over \mathrm{d} s}
      e_{ij}(\bm{u}(s))\,ds
      + \int_0^t {\mathrm{d} \over \mathrm{d} t}\widetilde\pi^P(t-s)\,p(s)\,ds
    """
    p_corr = convolve_field_scalar(corrs_pressure[vdp], pressures,
                                   iel, ts)
    p_corr += convolve_field_sym_tensor(corrs_rs, dstrains, vdp,
                                        dim, iel, ts)
    return p_corr


def compute_u_from_macro(strain, coor, iel, centre=None):
    r"""
    Macro-induced displacements.

    .. math::
      e_{ij}^x(\bm{u})\,(y_j - y_j^c)
    """
    n_nod, dim = coor.shape

    if centre is None:
        centre = nm.zeros((dim,), dtype=nm.float64)

    n_nod, dim = coor.shape
    um = nm.zeros((n_nod * dim,), dtype=nm.float64)
    for ir in range(dim):
        for ic in range(dim):
            ii = coor_to_sym(ir, ic, dim)
            um[ir::dim] += strain[iel, 0, ii, 0] * (coor[:, ic] - centre[ic])
    return um


def compute_p_from_macro(p_grad, coor, iel, centre=None, extdim=0):
    r"""
    Macro-induced pressure.

    .. math::
      \partial_j^x p\,(y_j - y_j^c)
    """
    n_nod, dim = coor.shape

    if centre is None:
        centre = nm.zeros((dim,), dtype=nm.float64)

    n_nod, dim = coor.shape
    pm = nm.zeros((n_nod,), dtype=nm.float64)
    for ic in range(dim + extdim):
        pm += p_grad[iel, 0, ic, 0] * (coor[:, ic] - centre[ic])
    return pm


def compute_micro_u(corrs, strain, vu, dim, out=None):
    r"""
    Micro displacements.

    .. math::
      \bm{u}^1 = \bm{\chi}^{ij}\, e_{ij}^x(\bm{u}^0)
    """

    if out is None:
        out = nm.zeros_like(corrs[vu+'_00'])

    for ir in range(dim):
        for ic in range(dim):
            ii = coor_to_sym(ir, ic, dim)
            out += corrs[vu+'_%d%d' % (ir, ic)] * strain[ii]
    return out


def compute_stress_strain_u(pb, integral, region, material, vu, data):
    var = pb.create_variables([vu])[vu]
    var.set_data(data)

    stress = pb.evaluate('ev_cauchy_stress.%s.%s(%s, %s)'
                         % (integral, region, material, vu), verbose=False,
                         mode='el_avg', **{vu: var})
    strain = pb.evaluate('ev_cauchy_strain.%s.%s(%s)'
                         % (integral, region, vu), verbose=False,
                         mode='el_avg', **{vu: var})

    return extend_cell_data(stress, pb.domain, region), \
           extend_cell_data(strain, pb.domain, region)


def add_stress_p(out, pb, integral, region, vp, data):
    var = pb.create_variables([vp])[vp]
    var.set_data(data)

    press0 = pb.evaluate('ev_integrate.%s.%s(%s)'
                         % (integral, region, vp), verbose=False,
                         mode='el_avg', **{vp: var})
    press = extend_cell_data(press0, pb.domain, region)

    dim = pb.domain.mesh.dim
    nn = out.shape[0]
    for ii in range(nn):
        for j in range(dim):
            out[ii, 0, j, 0] += press[ii, 0, 0, 0]


def compute_mac_stress_part(pb, integral, region, material, vu, mac_strain):
    avgmat = pb.evaluate('ev_integrate_mat.%s.%s(%s, %s)'
                         % (integral, region, material, vu), verbose=False,
                         mode='el_avg')

    return extend_cell_data(nm.dot(avgmat, mac_strain), pb.domain, region)


def recover_bones(problem, micro_problem, region, eps0,
                  ts, strain, dstrains, p_grad, pressures,
                  corrs_permeability, corrs_rs, corrs_time_rs,
                  corrs_pressure, corrs_time_pressure,
                  var_names, naming_scheme='step_iel'):
    r"""
    Notes
    -----
    - note that

      .. math::
        \widetilde{\pi}^P

      is in corrs_pressure -> from time correctors only 'u', 'dp' are needed.
    """

    dim = problem.domain.mesh.dim

    vu, vp, vn, vpp1, vppp1 = var_names
    vdp = 'd' + vp

    variables = micro_problem.create_variables()
    to_output = variables.create_output

    micro_u, micro_p = variables[vu], variables[vp]

    micro_coor = micro_u.field.get_coor()

    nodes_yc = micro_problem.domain.regions['Yc'].vertices

    join = os.path.join
    format = get_print_info(problem.domain.mesh.n_el, fill='0')[1]

    for ii, iel in enumerate(region.cells):
        print('ii: %d, iel: %d' % (ii, iel))

        pressure = pressures[-1][ii, 0, 0, 0]

        us = corrs_pressure[vu].data * pressure
        add_strain_rs(corrs_rs, strain, vu, dim, ii, out=us)

        ut = convolve_field_scalar(corrs_time_pressure[vu], pressures, ii, ts)
        ut += convolve_field_sym_tensor(corrs_time_rs, dstrains, vu,
                                        dim, ii, ts)
        u1 = us + ut

        u_mic = compute_u_from_macro(strain, micro_coor, ii) + u1

        ps = corrs_pressure[vp].data * pressure

        pt = convolve_field_scalar(corrs_time_pressure[vdp], pressures, ii, ts)
        pt += convolve_field_sym_tensor(corrs_time_rs, dstrains, vdp,
                                        dim, ii, ts)
        p_hat = ps + pt

        # \eta_k \partial_k^x p
        p1 = combine_scalar_grad(corrs_permeability, p_grad, vn, ii)

        p_hat_e = micro_p.field.extend_dofs(p_hat[:, None], fill_value=0.0)
        p_mic = compute_p_from_macro(p_grad, micro_coor, ii)[:, None] \
                + p_hat_e / eps0
        p_mic[nodes_yc] = p1[:, None]

        # (y_k + \eta_k) \partial_k^x p
        p_aux = combine_scalar_grad(corrs_permeability, p_grad, vn, ii,
                                    shift_coors=micro_coor[nodes_yc])

        meval = micro_problem.evaluate

        var_p = variables[vppp1]
        var_p.set_data(p_aux)
        dvel_m1 = meval('ev_diffusion_velocity.i1.Yc(m.K, %s)' % vppp1,
                        verbose=False, mode='el_avg', **{vppp1: var_p})

        var_p = variables[vpp1]
        var_p.set_data(p_hat)
        dvel_m2 = meval('ev_diffusion_velocity.i1.Ym(m.K, %s)' % vpp1,
                        verbose=False, mode='el_avg',
                        **{vpp1: var_p}) * eps0

        out = {}
        out.update(to_output(u_mic, var_info={vu: (True, vu)},
                             extend=True))
        out[vp] = Struct(name='output_data',
                         mode='vertex', data=p_mic,
                         var_name=vp, dofs=micro_p.dofs)

        aux = extend_cell_data(dvel_m1, micro_problem.domain, 'Yc')
        out['dvel_m1'] = Struct(name='output_data',
                                mode='cell', data=aux,
                                dofs=None)

        aux = extend_cell_data(dvel_m2, micro_problem.domain, 'Ym')
        out['dvel_m2'] = Struct(name='output_data',
                                mode='cell', data=aux,
                                dofs=None)

        suffix = get_output_suffix(iel, ts, naming_scheme, format,
                                   micro_problem.output_format)
        micro_name = micro_problem.get_output_name(extra=suffix)
        filename = join(problem.output_dir,
                        'recovered_' + os.path.basename(micro_name))

        micro_problem.save_state(filename, out=out, ts=ts)


def recover_paraflow(problem, micro_problem, region,
                     ts, strain, dstrains, pressures1, pressures2,
                     corrs_rs, corrs_time_rs,
                     corrs_alpha1, corrs_time_alpha1,
                     corrs_alpha2, corrs_time_alpha2,
                     var_names, naming_scheme='step_iel'):

    dim = problem.domain.mesh.dim

    vu, vp = var_names
    vdp = 'd' + vp

    micro_u = micro_problem.variables[vu]
    micro_coor = micro_u.field.get_coor()

    micro_p = micro_problem.variables[vp]

    nodes_y1 = micro_problem.domain.regions['Y1'].vertices
    nodes_y2 = micro_problem.domain.regions['Y2'].vertices

    to_output = micro_problem.variables.create_output

    join = os.path.join
    format = get_print_info(problem.domain.mesh.n_el, fill='0')[1]

    for ii, iel in enumerate(region.cells):
        print('ii: %d, iel: %d' % (ii, iel))

        p1, p2 = pressures1[-1][ii, 0, 0, 0], pressures2[-1][ii, 0, 0, 0]

        us = corrs_alpha1[vu].data * p1 + corrs_alpha2[vu].data * p2
        add_strain_rs(corrs_rs, strain, vu, dim, ii, out=us)

        ut = convolve_field_scalar(corrs_time_alpha1[vu], pressures1, ii, ts)
        ut += convolve_field_scalar(corrs_time_alpha2[vu], pressures2, ii, ts)
        ut += convolve_field_sym_tensor(corrs_time_rs, dstrains, vu,
                                        dim, ii, ts)

        u_corr = us + ut
        u_mic = compute_u_from_macro(strain, micro_coor, ii) + u_corr

        ps = corrs_alpha1[vp].data * p1 + corrs_alpha2[vp].data * p2

        pt = convolve_field_scalar(corrs_time_alpha1[vdp], pressures1, ii, ts)
        pt += convolve_field_scalar(corrs_time_alpha2[vdp], pressures2, ii, ts)
        pt += convolve_field_sym_tensor(corrs_time_rs, dstrains, vdp,
                                        dim, ii, ts)

        p_corr = ps + pt

        p_mic = micro_p.field.extend_dofs(p_corr[:, nm.newaxis])
        p_mic[nodes_y1] = p1
        p_mic[nodes_y2] = p2

        out = {}
        out.update(to_output(u_mic, var_info={vu: (True, vu)}, extend=True))
        out[vp] = Struct(name='output_data',
                         mode='vertex', data=p_mic,
                         var_name=vp, dofs=micro_p.dofs)

        suffix = get_output_suffix(iel, ts, naming_scheme, format,
                                   micro_problem.output_format)
        micro_name = micro_problem.get_output_name(extra=suffix)
        filename = join(problem.output_dir, 'recovered_' + micro_name)

        micro_problem.save_state(filename, out=out, ts=ts)


def save_recovery_region(mac_pb, rname, filename=None):
    filename = get_default(filename, os.path.join(mac_pb.output_dir,
                                                  'recovery_region.vtk'))

    region = mac_pb.domain.regions[rname]

    # Save recovery region characteristic function.
    out = {}
    mask = region.get_charfun(by_cell=False, val_by_id=False)
    out['vmask'] = Struct(name='output_data',
                          mode='vertex', data=mask[:, nm.newaxis],
                          dofs=None)
    mask = region.get_charfun(by_cell=True, val_by_id=False)
    out['cmask'] = Struct(name='output_data',
                          mode='cell',
                          data=mask[:, nm.newaxis, nm.newaxis, nm.newaxis],
                          dofs=None)

    mac_pb.save_state(filename, out=out)


_recovery_global_dict = {}


def destroy_pool():
    if 'pool' in _recovery_global_dict:
        _recovery_global_dict['pool'].close()


atexit.register(destroy_pool)


def _recovery_hook(args):
    idx, label, local_macro, verbose = args
    pb, corrs, recovery_hook = _recovery_global_dict['hook_args']

    output.set_output(quiet=False)
    output.level = label[1]
    output(label[0])
    output.set_output(quiet=not(verbose))

    if isinstance(corrs, list):
        corrs = corrs[idx]

    return recovery_hook(pb, corrs, local_macro)


def get_recovery_points(region, eps0):
    rcoors = region.domain.mesh.coors[region.get_entities(0), :]
    rcmin = nm.min(rcoors, axis=0)
    rcmax = nm.max(rcoors, axis=0)
    nn = nm.round((rcmax - rcmin) / eps0)
    if nm.prod(nn) == 0:
        raise ValueError(
            'incompatible recovery region and microstructure size!')

    cs = []
    for ii, n in enumerate(nn):
        cs.append(nm.arange(n) * eps0 + rcmin[ii])

    ccoors = nm.empty((int(nm.prod(nn)), nn.shape[0]), dtype=nm.float64)
    for ii, icoor in enumerate(nm.meshgrid(*cs, indexing='ij')):
        ccoors[:, ii] = icoor.flatten()

    return ccoors


def recover_micro_hook(micro_filename, region, macro, eps0,
                       region_mode='el_centers', eval_mode='constant',
                       eval_vars=None, corrs=None, recovery_file_tag='',
                       define_args=None, output_dir=None, verbose=False):
    """
    Parameters
    ----------
    micro_filename : str
        The definition file of the microproblem.
    region : Region or array
        The macroscopic region to be recovered. If array, the centers of
        microscopic RVEs (Representative Volume Element) are assumed to be
        stored in it. If Region, the RVE centers are computed according to
        `region_mode`, see below.
    macro : dict of arrays or tuples
        Either macroscopic values (if array) or the tuple
        (mode, eval_var, nodal_values) is expected. The tuple is used to
        evaluate the macroscopic values in given points of RVEs (see
        'eval_mode`). `mode` can be 'val', 'grad', 'div', or 'cauchy_strain'.
    eps0 : float
        The size of the microstructures (RVE).
    region_mode : {'el_centers', 'tiled'}
        If 'el_centers', the RVE centers are identical to the element centers
        of the macroscopic FE mesh. If 'tiled', the recovered region is tiled
        by rescaled RVEs.
    eval_mode : {'constant', 'continuous'}
        If 'constant', the macroscopic fields are evaluated only at the RVE
        centers. If 'continuous', the fields are evaluated at all points of
        the RVE mesh.
    eval_vars : list of variables
        The list of variables use to evaluate the macroscopic fields.
    corrs : dict of CorrSolution
        The correctors for recovery.
    recovery_file_tag : str
        The tag which is appended to the output file.
    define_args : dict
        The define arguments for the microscopic problem.
    output_dir : str
        The output directory.
    verbose : bool
        The verbose terminal output.
    """
    import sfepy.base.multiproc_proc as multi

    if 'micro_problem' not in _recovery_global_dict:
        # Create a micro-problem instance.
        required, other = get_standard_keywords()
        required.remove('equations')
        conf = ProblemConf.from_file(micro_filename, required, other,
                                     verbose=False, define_args=define_args)
        if output_dir is not None:
            conf.options.output_dir = output_dir

        recovery_hook = conf.options.get('recovery_hook', None)
        pb = None
        if recovery_hook is not None:
            recovery_hook = conf.get_function(recovery_hook)

            if corrs is None:
                # Coefficients and correctors
                coefs_filename = conf.options.get('coefs_filename', 'coefs')
                coefs_filename = op.join(conf.options.get('output_dir', '.'),
                                         coefs_filename + '.h5')
                coefs = Coefficients.from_file_hdf5(coefs_filename)
                corrs = \
                    get_correctors_from_file_hdf5(dump_names=coefs.save_names)

            pb = Problem.from_conf(conf, init_equations=False,
                                   init_solvers=False)

        _recovery_global_dict['micro_problem'] = pb, corrs, recovery_hook
    else:
        pb, corrs, recovery_hook = _recovery_global_dict['micro_problem']

    is_multiproc = pb.conf.options.get('multiprocessing', True)\
        and multi.use_multiprocessing

    if recovery_hook is not None:
        if eval_mode not in ['constant', 'continuous']:
            raise ValueError(f"Evaluation mode '{eval_mode}' not implemented!")

        if isinstance(region, Region):
            if region_mode == 'el_centers':
                dim = region.domain.mesh.dim
                ccoors = region.cmesh.get_centroids(dim)
                ccoors = ccoors[region.get_cells(), :]
                recovery_ids = region.get_entities(-1)
                recovery_ids_s = [f'(element {k})' for k in recovery_ids]
            elif region_mode == 'tiled':
                ccoors = get_recovery_points(region, eps0)
                recovery_ids = nm.arange(ccoors.shape[0])
                recovery_ids_s = [f'(point {k})' for k in recovery_ids]
            else:
                raise ValueError(
                    f"Region mode '{region_mode}' not implemented!")

        else:
            ccoors = region
            recovery_ids = nm.arange(ccoors.shape[0])
            recovery_ids_s = [f'(point {k})' for k in recovery_ids]

        nrec = ccoors.shape[0]

        output(f'recovering {nrec} microsctructures...')
        timer = Timer(start=True)
        output_level, output_fun = output.level, output.output_function

        _recovery_global_dict['hook_args'] = pb, corrs, recovery_hook

        if is_multiproc:
            multi_local_macro = []
            num_workers = nm.min([multi.cpu_count(), nrec])
            _recovery_global_dict['pool'] = multi.Pool(processes=num_workers)
        else:
            outs = []

        # Recover micro. region
        mesh = pb.domain.mesh
        bbox = mesh.get_bounding_box()
        mic_coors = (mesh.coors - 0.5 * (bbox[1, :] + bbox[0, :])) * eps0
        coors, conn, ndoffset, rec_ids = [], [], 0, []

        if eval_vars is not None and eval_mode == 'constant':
            new_macro = {}
            for k, v in macro.items():
                if isinstance(v, tuple):
                    emode, evar, val = v
                    efield = eval_vars[evar].field
                    new_macro[k] = efield.evaluate_at(ccoors, val, mode=emode)
                else:
                    new_macro[k] = v

            macro = new_macro

        for ii in range(nrec):
            local_coors = mic_coors.copy()
            local_coors[:, :ccoors.shape[1]] += ccoors[ii]
            coors.append(local_coors)
            conn.append(mesh.get_conn(mesh.descs[0]) + ndoffset)
            ndoffset += mesh.n_nod
            rec_ids.append(nm.ones((mesh.n_el,)) * recovery_ids[ii])

            label = (f'micro: {ii + 1}/{nrec} {recovery_ids_s[ii]}',
                     output_level)

            local_macro = {'eps0': eps0}

            if eval_vars is not None and eval_mode == 'continuous':
                for k, v in macro.items():
                    if isinstance(v, tuple):
                        emode, evar, val = v
                        efield = eval_vars[evar].field

                        # # Inside recovery region?
                        # from sfepy.base.base import debug; debug()
                        # v = nm.ones((efield.region.entities[0].shape[0], 1))
                        # v[evfield.vertex_remap[region.entities[0]]] = 0
                        # no = nm.sum(v)
                        # aux = efield.evaluate_at(local_coors, v)
                        # if no > 0 and (nm.sum(aux) / no) > 1e-3:
                        #     print('not inside!')
                        #     continue

                        local_macro[k] = efield.evaluate_at(local_coors, val,
                                                            mode=emode)
                    else:
                        local_macro[k] = v
            else:
                for k, v in macro.items():
                    if k.startswith('_'):
                        local_macro[k] = v
                    else:
                        local_macro[k] = v[ii, ...]

            if is_multiproc:
                multi_local_macro.append((ii, label, local_macro, verbose))
            else:
                outs.append(_recovery_hook((ii, label, local_macro, verbose)))

        rec_ids = nm.hstack(rec_ids).reshape((-1, 1, 1, 1))

        if is_multiproc:
            pool = _recovery_global_dict['pool']
            outs = list(pool.map(_recovery_hook, multi_local_macro))

        output.level, output.output_function = output_level, output_fun

        # Collect output data - by region
        outregs_data = {}
        outregs_info = {None: ('ALL', nm.arange(pb.domain.mesh.n_el))}
        vn_reg = {None: None}
        for k, v in outs[0].items():
            if hasattr(v, 'region_name'):
                rn = v.region_name
                if rn not in outregs_info:
                    reg = pb.domain.regions[rn]
                    outregs_info[rn] = (rn, reg.get_entities(-1))
            else:
                vn = getattr(v, 'var_name', None)
                if vn not in vn_reg:
                    reg = pb.create_variables(vn)[vn].field.region
                    rn = reg.name
                    vn_reg[vn] = rn
                    if rn not in outregs_info:
                        outregs_info[rn] = (rn, reg.get_entities(-1))
                else:
                    rn = vn_reg[vn]

            if rn in outregs_data:
                outregs_data[rn].append(k)
            else:
                outregs_data[rn] = [k]

        nrve = len(coors)
        coors = nm.vstack(coors)
        ngroups = nm.tile(mesh.cmesh.vertex_groups.squeeze(), (nrve,))
        conn = nm.vstack(conn)
        cgroups = nm.tile(mesh.cmesh.cell_groups.squeeze(), (nrve,))
        # Get region mesh and data
        output_dir = pb.conf.options.get('output_dir', '.')
        for rn in outregs_data.keys():
            rlabel, cidxs = outregs_info[rn]
            gcidxs = nm.hstack([cidxs + mesh.n_el * ii for ii in range(nrve)])
            rconn = conn[gcidxs]
            remap = -nm.ones((coors.shape[0],), dtype=nm.int32)
            remap[rconn] = 1
            vidxs = nm.where(remap > 0)[0]
            remap[vidxs] = nm.arange(len(vidxs))
            rconn = remap[rconn]
            rcoors = coors[vidxs, :]

            out = {}
            for ivar in outregs_data[rn]:
                data = [outs[ii][ivar].data for ii in range(nrve)]
                out[ivar] = Struct(name='output_data',
                                   mode=outs[0][ivar].mode,
                                   data=nm.vstack(data))

            out['recovery_id'] = Struct(name='output_data',
                                        mode='cell',
                                        data=rec_ids[gcidxs, ...])
            micro_name = pb.get_output_name(extra='recovered_%s%s'
                                            % (rlabel, recovery_file_tag))
            filename = op.join(output_dir, op.basename(micro_name))
            mesh_out = Mesh.from_data('recovery_%s' % rlabel, rcoors,
                                      ngroups[vidxs], [rconn],
                                      [cgroups[gcidxs]], [mesh.descs[0]])
            mesh_out.write(filename, io='auto', out=out)
            output(f'output saved to "{filename}"')

        output('...done in %.2f s' % timer.stop())

import numpy as nm
import numpy.linalg as nla

from sfepy.fem.mappings import get_physical_qps

def eval_ion_ion_energy(centres, charges):
    """
    Compute the ion-ion nergy.
    """
    val = 0.0
    for ir, rcentre in enumerate(centres):
        for ic, ccentre in enumerate(centres):
            if ir == ic: continue

            val += charges[ir] * charges[ic] / nla.norm(rcentre - ccentre)

    val *= 0.5

    return val

def eval_non_local_interaction(problem, region_name, var_name,
                               integral_name, f1, f2, kernel_function,
                               pbar=None):
    """
    Single element group only!
    """
    var = problem.get_variables()[var_name]

    region = problem.domain.regions[region_name]

    integral = problem.integrals[integral_name]

    qps = get_physical_qps(region, integral)
    igs = qps.values.keys()

    ig = igs[0]
    qp = qps.values[ig]

    n_el, n_qp = f1.shape[:2]

    vg, _ = var.get_mapping(ig, region, integral, 'volume')

    # Weighted jacobian.
    det = vg.variable(1)

    shape = (n_el * n_qp, 1, 1)

    val1 = f1 * det
    val2 = f2 * det

    val1.shape = shape
    val2.shape = shape

    coef = nm.zeros(shape, dtype=val1.dtype)

    if pbar is not None:
        pbar.init(qp.shape[0])
        pbar.update(0)

    for ii, coor in enumerate(qp):
        ## tt = time.clock()
        kernel = kernel_function(coor, qp)
        kernel.shape = val2.shape
        ## print'aa',   time.clock() - tt

        ## tt = time.clock()
        coef[ii] = (kernel * val2).sum()
        ## print 'bb', time.clock() - tt

        if pbar is not None:
            pbar.update(ii)

    val = (val1 * coef).sum()

    if pbar is not None:
        print

    return val

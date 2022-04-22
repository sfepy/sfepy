from __future__ import absolute_import
import numpy as nm
import scipy.sparse as sp
from scipy.sparse.linalg import aslinearoperator
from scipy.linalg import eigvals_banded

from sfepy.base.base import get_default, output
from sfepy.linalg import infinity_norm
from six.moves import range

def sym_tri_eigen(diags, select_indices=None):
    """
    Compute eigenvalues of a symmetric tridiagonal matrix using
    `scipy.linalg.eigvals_banded()`.
    """
    if select_indices is not None:
        n = diags.shape[1]
        select_indices = nm.minimum(select_indices, n)
        eigs = eigvals_banded(diags, lower=True, select='i',
                              select_range=select_indices)

    else:
        eigs = eigvals_banded(diags, lower=True, select='a')

    return eigs

def cg_eigs(mtx, rhs=None, precond=None, i_max=None, eps_r=1e-10,
            shift=None, select_indices=None, verbose=False, report_step=10):
    r"""
    Make several iterations of the conjugate gradients and estimate so
    the eigenvalues of a (sparse SPD) matrix (Lanczos algorithm).

    Parameters
    ----------
    mtx : spmatrix or array
        The sparse matrix :math:`A`.
    precond : spmatrix or array, optional
        The preconditioner matrix. Any object that can be multiplied by
        vector can be passed.
    i_max : int
        The maximum number of the Lanczos algorithm iterations.
    eps_r : float
        The relative stopping tolerance.
    shift : float, optional
        Eigenvalue shift for non-SPD matrices. If negative, the shift is
        computed as :math:`|shift| ||A||_{\infty}`.
    select_indices : (min, max), optional
        If given, computed only the eigenvalues with indices `min <= i <= max`.
    verbose : bool
        Verbosity control.
    report_step : int
        If `verbose` is True, report in every `report_step`-th step.

    Returns
    -------
    vec : array
        The approximate solution to the linear system.
    n_it : int
        The number of CG iterations used.
    norm_rs : array
        Convergence history of residual norms.
    eigs : array
        The approximate eigenvalues sorted in ascending order.
    """
    n_row = mtx.shape[0]
    norm = nm.linalg.norm

    rhs = get_default(rhs, nm.random.rand(n_row))
    i_max = get_default(i_max, min(n_row, 100))

    if shift is not None:
        if shift < 0.0:
            mtx_norm = infinity_norm(mtx)
            output('matrix max norm:', mtx_norm, verbose=verbose)
            shift = abs(shift) * mtx_norm

        output('eigenvalue shift:', shift, verbose=verbose)
        mtx = mtx + shift * sp.eye(n_row, n_row, dtype=mtx.dtype)

    mtx = aslinearoperator(mtx)

    lambda_max = 0
    lambda_min = 0
    econd = 1.0

    x0 = nm.zeros_like(rhs)
    r = r0 = rhs

    # The diagonals (0, 1) in two rows.
    diags = nm.empty((2, i_max + 1), dtype=mtx.dtype)
    diags[0, 0] = 0

    if precond is None:
        z0 = r0

    else:
        z0 = precond * r0

    x = x0
    z = z0

    p = nm.zeros_like(rhs)
    beta = 0.0

    rho0 = nm.dot(z0, r0)
    norm_r0 = norm(r0);

    if verbose:
        output('%5d lambda: %e %e cond: %e |R|: %e\n' % (0, 0, 0, 0, norm_r0))

    norm_rs = [norm_r0]

    for ii in range(i_max):
        p = z + beta * p
        q = mtx * p

        alpha = rho0 / nm.dot(p, q)
        if (not nm.isfinite(alpha)) or abs(alpha) < 1e-16:
            output('precision limit reached!')
            ii -= 1
            break

        x = x + alpha * p
        r = r - alpha * q

        if precond is None:
            z = r

        else:
            z = precond * r

        norm_r = norm(r)
        norm_rs.append(norm_r)

        rho1 = nm.dot(z, r)
        beta = rho1 / rho0

        # Lanczos
        diags[0, ii] += 1.0 / alpha
        diags[1, ii] = - nm.sqrt(beta) / alpha
        diags[0, ii+1] = beta / alpha

        if verbose and (ii > 0):
            if (ii % report_step) == 0:
                eigs = sym_tri_eigen(diags[:, :ii+1],
                                     select_indices=select_indices)
                if select_indices is None:
                    lambda_min, lambda_max = eigs[0], eigs[-1]
                    econd = lambda_max / lambda_min
                    output('%5d lambda: %e %e cond: %e |R|: %e\n'
                           % (ii, lambda_min, lambda_max, econd, norm_r))

                else:
                    output('%5d |R|: %e\n'
                           % (ii, norm_r))

        if (norm_r / norm_r0) < eps_r:
            output('converged on %d iters, |Ri|/|R0|: %e, econd: %e\n'
                   % (ii, norm_r / norm_r0, econd), verbose=verbose)
            break

        rho0 = rho1

    eigs = sym_tri_eigen(diags[:, :ii+1], select_indices=select_indices)
    if verbose and (select_indices is None):
        lambda_min, lambda_max = eigs[0], eigs[-1]
        econd = lambda_max / lambda_min
        output('min: %e  max: %e  cond: %e\n'
               % (lambda_min, lambda_max, econd))

    if shift is not None:
        eigs -= shift

    return x, ii, nm.array(norm_rs), eigs

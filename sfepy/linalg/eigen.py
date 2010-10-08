import numpy as nm
import scipy.sparse as sps

from sfepy.base.base import get_default, output

def sym_tri_eigen(diags):
    """
    Compute eigenvalues of a symmetric tridiagonal matrix. Naive
    implementation!
    """
    n_row = diags.shape[1]
    shape = (n_row, n_row)

    aux = nm.empty((3, n_row + 1), dtype=diags.dtype)
    aux[:2, :-1] = diags
    aux[2, 1:] = diags[1]
    mtx = sps.dia_matrix((aux, (0, -1, 1)), shape=shape)

    eigs, _ = nm.linalg.eig(mtx.toarray())
    eigs.sort()

    return eigs

def cg_eigs(mtx, rhs=None, precond=None, i_max=None, eps_a=1e-10,
            verbose=False, report_step=10):
    """
    Make several iterations of the conjugate gradients and estimate so
    the eigenvalues of a (sparse) SPD matrix (Lanczos algorithm).

    Returns
    -------
    vec : array
        The approximate solution to the linear system.
    n_it : int
        The number of CG iterations used.
    norm_rs : array
        Convergence history of residual norms
    eigs : array
        The approximate eigenvalues sorted in ascending order.
    """
    n_row = mtx.shape[0]
    norm = nm.linalg.norm

    rhs = get_default(rhs, nm.random.rand(n_row))
    i_max = get_default(i_max, min(n_row, 100))

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

    for ii in xrange(i_max):
        p = z + beta * p
        q = mtx * p

        alpha = rho0 / nm.dot(p, q)

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

        if verbose:
            if (ii % report_step) == 1:
                eigs = sym_tri_eigen(diags[:, :ii+1])
                lambda_min, lambda_max = eigs[0], eigs[-1]
                econd = lambda_max / lambda_min
                output('%5d lambda: %e %e cond: %e |R|: %e\n'
                       % (ii, lambda_min, lambda_max, econd, norm_r))

        if (norm_r / norm_r0) < eps_a:
            if verbose:
                output('converged on %d iters, |Ri|/|R0|: %e, econd: %e\n'
                       % (ii, norm_r / norm_r0, econd))
            break

        rho0 = rho1

    eigs = sym_tri_eigen(diags[:, :ii+1])
    lambda_min, lambda_max = eigs[0], eigs[-1]
    econd = lambda_max / lambda_min
    if verbose:
        output('min: %e  max: %e  cond: %e\n'
               % (lambda_min, lambda_max, econd))

    return x, ii, nm.array(norm_rs), eigs

def arpack_eigs(mtx, nev=1, which='SM'):
    """
    Calculate several eigenvalues and corresponding eigenvectors of a
    matrix using ARPACK from SciPy. The eigenvalues are sorted in
    ascending order.
    """
    from scipy.sparse.linalg.interface import aslinearoperator
    from scipy.sparse.linalg.eigen.arpack import speigs

    matvec = aslinearoperator(mtx).matvec
    eigs, vecs = speigs.ARPACK_eigs(matvec, mtx.shape[0], nev=nev, which=which)

    ii = nm.argsort(eigs)
    eigs = eigs[ii]
    vecs = vecs[:,ii]

    return eigs, vecs

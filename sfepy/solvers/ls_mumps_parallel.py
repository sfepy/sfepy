import numpy as nm
from mpi4py import MPI
import ls_mumps as mumps
from tempfile import gettempdir
import os.path as op


def tmpfile(fname):
    return op.join(gettempdir(), fname)


def mumps_parallel_solve():
    comm = MPI.COMM_WORLD
    flags = nm.memmap('vals_falgs.array', dtype='int32', mode='r')

    mumps_ls = mumps.MumpsSolver(system={0: 'real', 1: 'complex'}[flags[1]])
    if comm.rank == 0:
        dtype = {0: 'float64', 1: 'complex128'}[flags[1]]

        vals_mtx = nm.memmap(tmpfile('vals_mtx.array'), dtype=dtype, mode='r')
        vals_rhs = nm.memmap(tmpfile('vals_rhs.array'), dtype=dtype, mode='r')
        idxs = nm.memmap(tmpfile('idxs.array'), dtype='int32',
                         mode='r').reshape((2, vals_mtx.shape[0]))

        n = flags[0]
        mumps_ls.set_rcd_centralized(idxs[0, :], idxs[1, :], vals_mtx, n)
        x = vals_rhs.copy()
        mumps_ls.set_rhs(x)

    mumps_ls(6)  # analyse, factorize, solve

    if comm.rank == 0:
        vals_x = nm.memmap(tmpfile('vals_x.array'), dtype=dtype, mode='w+',
                           shape=x.shape)
        vals_x[:] = x

    del(mumps_ls)


comm = MPI.Comm.Get_parent()
mumps_parallel_solve()
comm.Disconnect()

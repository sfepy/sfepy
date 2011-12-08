#!/usr/bin/env python
"""
PETSc solver worker process.
"""
import time
import sys, petsc4py
petsc4py.init(sys.argv)

from petsc4py import PETSc

def solve():
    opts = PETSc.Options()

    mat_filename = opts.getString('-mtx', 'mtx.dat')
    rhs_filename = opts.getString('-rhs', 'rhs.dat')
    sol0_filename = opts.getString('-sol0', '')
    sol_filename = opts.getString('-sol', 'sol.dat')
    status_filename = opts.getString('-status', 'status.txt')

    view_mtx = PETSc.Viewer().createBinary(mat_filename, mode='r')
    view_rhs = PETSc.Viewer().createBinary(rhs_filename, mode='r')

    mtx = PETSc.Mat().load(view_mtx)
    rhs = PETSc.Vec().load(view_rhs)

    if not sol0_filename:
        sol = rhs.duplicate()

    else:
        view_sol0 = PETSc.Viewer().createBinary(sol0_filename, mode='r')
        sol = PETSc.Mat().load(view_sol0)

    ksp = PETSc.KSP().create()
    ksp.setOperators(mtx)
    ksp.setFromOptions()

    tt = time.clock()
    ksp.solve(rhs, sol)
    elapsed = time.clock() - tt

    view_sol = PETSc.Viewer().createBinary(sol_filename, mode='w')
    sol.view(view_sol)

    fd = open(status_filename, 'w')
    fd.write('%d %.2f' % (ksp.reason, elapsed))
    fd.close()

if __name__ == '__main__':
    solve()

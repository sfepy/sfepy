#! /usr/bin/python
from pysparse import spmatrix, jdsym, itsolvers, precon
from numpy import array

def assemble(fname,topo,mapping):
    nnodes=len(mapping)
    f=file(fname)
    A=spmatrix.ll_mat(nnodes,nnodes)
    M=spmatrix.ll_mat(nnodes,nnodes)
    for el in topo:
        n=len(el)
        matK=readmat(f,n)
        matV=readmat(f,n)
        matM=readmat(f,n)
        el=[mapping[x] for x in el]
        A.update_add_mask_sym(matK,array(el),ones(n))
        A.update_add_mask_sym(matV,array(el),ones(n))
        M.update_add_mask_sym(matM,array(el),ones(n))
    return A,M

def convert_mat(self):
    A = spmatrix.ll_mat(*self.shape)
    for ii in xrange(self.size):
        i, j = self.rowcol(ii)
        val = self.getdata(ii)
        A.update_add_at([val], [i], [j])
    return A

def read_mat(fname):
    f = open(fname)
    nnodes = int(f.readline().split(" ")[0])
    A = spmatrix.ll_mat(nnodes, nnodes)
    nlines = int(f.readline().split(" ")[0])
    for l in range(nlines):
        i, j, val = f.readline().split(" ")
        i = int(i); j = int(j); val = float(val)
        A.update_add_at([val], [i], [j])
    return A

def savevec(vec,fname):
    f=file(fname,"w")
    for n in vec:
        f.write("%f "%n)
    f.write("\n")

def solve(A, B, nE = 10):
    print "loading..."
    A = convert_mat(A)
    M = convert_mat(B)
    print "solving..."
    tau=0.0
    #tau=-10
    Atau=A.copy()
    Atau.shift(-tau,M)
    K=precon.jacobi(Atau)
    #nE=100
    #K=None
    A=A.to_sss();M=M.to_sss();
    kconv, lmbd, Q, it, it_in = jdsym.jdsym(A, M, K, nE, tau, 1e-5, 150, 
            itsolvers.qmrs, clvl=1, strategy=1)

    print "number of converged eigenvalues:",kconv
    print "done"
    return lmbd, Q
    #print lmbd
    #for i in range(kconv):
    #    savevec(Q[:,i],"tmp/sol-%d.dat"%i)


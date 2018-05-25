import ctypes
try:
    from mpi4py import MPI
    use_mpi = True
except ImportError:
    use_mpi = False

MUMPS_VERSION_MAX_LEN = 14

ctypes = ctypes
c_pointer = ctypes.POINTER
c_int = ctypes.c_int
c_double = ctypes.c_double

mumps_int = c_int
mumps_pint = c_pointer(mumps_int)
mumps_real = c_double
mumps_preal = c_pointer(mumps_real)
mumps_complex = c_double
mumps_pcomplex = c_pointer(mumps_complex)

mumps_libs = {}


def load_library(libname):
    """Load shared library in a system dependent way."""
    import sys

    if sys.platform.startswith('win'):  # Windows system
        from ctypes.util import find_library

        lib_fname = find_library(libname)

    else:  # Linux system
        lib_fname = 'lib' + libname + '.so'

    lib = ctypes.cdll.LoadLibrary(lib_fname)

    return lib


def load_mumps_libraries():
    mumps_libs['dmumps'] = load_library('dmumps').dmumps_c
    mumps_libs['zmumps'] = load_library('zmumps').zmumps_c


class mumps_struc_c(ctypes.Structure):
    # MUMPS 4.10.0
    _fields_ = [
        ('sym', mumps_int),
        ('par', mumps_int),
        ('job', mumps_int),
        ('comm_fortran', mumps_int),
        ('icntl', mumps_int * 40),
        ('cntl', mumps_real * 15),
        ('n', mumps_int),
        #
        ('nz_alloc', mumps_int),
        # /* Assembled entry */
        ('nz', mumps_int),
        ('irn', mumps_pint),
        ('jcn', mumps_pint),
        ('a', mumps_pcomplex),
        # /* Distributed entry */
        ('nz_loc', mumps_int),
        ('irn_loc', mumps_pint),
        ('jcn_loc', mumps_pint),
        ('a_loc', mumps_pcomplex),
        # /* Element entry */
        ('nelt', mumps_int),
        ('eltptr', mumps_pint),
        ('eltvar', mumps_pint),
        ('a_elt', mumps_pcomplex),
        # /* Ordering, if given by user */
        ('perm_in', mumps_pint),
        # /* Orderings returned to user */
        ('sym_perm', mumps_pint),
        ('uns_perm', mumps_pint),
        # /* Scaling (input only in this version) */
        ('colsca', mumps_preal),
        ('rowsca', mumps_preal),
        # /* RHS, solution, ouptput data and statistics */
        ('rhs', mumps_pcomplex),
        ('redrhs', mumps_pcomplex),
        ('rhs_sparse', mumps_pcomplex),
        ('sol_loc', mumps_pcomplex),
        ('irhs_sparse', mumps_pint),
        ('irhs_ptr', mumps_pint),
        ('isol_loc', mumps_pint),
        ('nrhs', mumps_int),
        ('lrhs', mumps_int),
        ('lredrhs', mumps_int),
        ('nz_rhs', mumps_int),
        ('isol_loc', mumps_int),
        ('schur_mloc', mumps_int),
        ('schur_nloc', mumps_int),
        ('schur_lld', mumps_int),
        ('mblock', mumps_int),
        ('nblock', mumps_int),
        ('nprow', mumps_int),
        ('npcol', mumps_int),
        ('info', mumps_int * 40),
        ('infog', mumps_int * 40),
        ('rinfo', mumps_real * 40),
        ('rinfog', mumps_real * 40),
        # /* Null space */
        ('deficiency', mumps_int),
        ('pivnul_list', mumps_pint),
        ('mapping', mumps_pint),
        # /* Schur */
        ('size_schur', mumps_int),
        ('listvar_schur', mumps_pint),
        ('schur', mumps_pcomplex),
        # /* Internal parameters */
        ('instance_number', mumps_int),
        ('wk_user', mumps_pcomplex),
        # /* Version number:
        #  length=14 in FORTRAN + 1 for final \0 + 1 for alignment */
        ('version_number', ctypes.c_char * (MUMPS_VERSION_MAX_LEN + 1 + 1)),
        # /* For out-of-core */
        ('ooc_tmpdir', ctypes.c_char * 256),
        ('ooc_prefix', ctypes.c_char * 64),
        # /* To save the matrix in matrix market format */
        ('write_problem', ctypes.c_char * 256),
        ('lwk_user', mumps_int)]


class MumpsSolver(object):
    """MUMPS object."""

    def __init__(self, sym=0, mpi_comm=None, system='real', silent=True):
        """
        Init MUMUPS solver.

        Parameters
        ----------
        sym : int
            1 = symmetric system
        mpi_comm : MPI Communicator or None
            If None, use MPI.COMM_WORLD
        system : one of 'real', 'comples'
            Use real or complex linear solver.
        silent : bool
            If True, no MUMPS error, warning, and diagnostic messages.
        """
        self.struct = None

        if not use_mpi:
            raise AttributeError('No mpi4py found, required by MUMPS solver!')

        if len(mumps_libs) == 0:
            load_mumps_libraries()

        if system == 'real':
            self._mumps_c = mumps_libs['dmumps']
        elif system == 'complex':
            self._mumps_c = mumps_libs['zmumps']

        self._mumps_c.argtypes = [c_pointer(mumps_struc_c)]
        self._mumps_c.restype = None

        self.struct = mumps_struc_c()
        self.struct.par = 1
        self.struct.sym = sym
        self.struct.n = 0
        self._mumps_call(job=-1)  # init
        self._data = {}

        self.mpi_comm = MPI.COMM_WORLD if mpi_comm is None else mpi_comm
        self.struct.comm_fortran = self.mpi_comm.py2f()
        self.rank = self.mpi_comm.rank

        # be silent
        if silent:
            self.set_silent()

    def set_silent(self):
        self.struct.icntl[0] = -1  # suppress error messages
        self.struct.icntl[1] = -1  # suppress diagnostic messages
        self.struct.icntl[2] = -1  # suppress global info
        self.struct.icntl[3] = 0

    def set_verbose(self):
        self.struct.icntl[0] = 6  # error messages
        self.struct.icntl[1] = 0  # diagnostic messages
        self.struct.icntl[2] = 6  # global info
        self.struct.icntl[3] = 2

    def __del__(self):
        """Finish MUMPS."""
        if self.struct is not None:
            self._mumps_call(job=-2)  # done

        self.struct = None

    def set_A_centralized(self, A):
        """
        Set the sparse matrix.

        Parameters
        ----------
        A : scipy sparse martix
            The sparse matrix.
        """
        # if not(self.rank == 0):
        #     return

        assert A.shape[0] == A.shape[1]

        A = A.tocoo()
        rr = A.row + 1
        cc = A.col + 1
        self.set_rcd_centralized(rr, cc, A.data, A.shape[0])

    def set_rcd_centralized(self, ir, ic, data, n):
        """
        Set the matrix by row and column indicies and data vector.
        The matrix shape is determined by the maximal values of
        row and column indicies.

        Parameters
        ----------
        ir : array
            The row idicies.
        ic : array
            The column idicies.
        data : array
            The matrix entries.
        n : int
            The matrix dimension.
        """
        # if not(self.rank == 0):
        #     return

        assert ir.shape[0] == ic.shape[0] == data.shape[0]

        self._data.update(ir=ir, ic=ic, data=data)
        self.struct.n = n
        self.struct.nz = ir.shape[0]
        self.struct.irn = ir.ctypes.data_as(mumps_pint)
        self.struct.jcn = ic.ctypes.data_as(mumps_pint)
        self.struct.a = data.ctypes.data_as(mumps_pcomplex)

    def set_b(self, b):
        """Set the right hand side of the linear system."""
        self._data.update(b=b)
        self.struct.rhs = b.ctypes.data_as(mumps_pcomplex)

    def __call__(self, job):
        """Set the job and call MUMPS."""
        self._mumps_call(job)

    def _mumps_call(self, job):
        """Set the job and call MUMPS.

        Jobs:
        -----
        1: analyse
        2: factorize
        3: solve
        4: analyse, factorize
        5: factorize, solve
        6: analyse, factorize, solve
        """
        self.struct.job = job
        self._mumps_c(ctypes.byref(self.struct))
        if self.struct.infog[0] < 0:
            raise RuntimeError("MUMPS error: %d" % self.struct.infog[0])

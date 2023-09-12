import ctypes
import re
import numpy as nm
try:
    from mpi4py import MPI
    use_mpi = True
except ImportError:
    use_mpi = False

AUX_LENGTH = 16 * 1024

c_pointer = ctypes.POINTER

mumps_int = ctypes.c_int
mumps_pint = c_pointer(mumps_int)
mumps_int8 = ctypes.c_uint64
mumps_real = ctypes.c_double
mumps_preal = c_pointer(mumps_real)
mumps_complex = ctypes.c_double
mumps_pcomplex = c_pointer(mumps_complex)

mumps_libs = {}

def dec(val, encoding='utf-8'):
    """
    Decode given bytes using the specified encoding.
    """
    import sys
    if isinstance(val, bytes) and sys.version_info > (3, 0):
        val = val.decode(encoding)
    return val


def load_library(libname):
    """Load shared library in a system dependent way."""
    import sys

    if sys.platform.startswith('win'):  # Windows system
        from ctypes.util import find_library

        lib_fname = find_library(libname)
        if lib_fname is None:
            lib_fname = find_library('lib' + libname)

    else:  # Linux system
        lib_fname = 'lib' + libname + '.so'

    lib = ctypes.cdll.LoadLibrary(lib_fname)

    return lib


def load_mumps_libraries():
    mumps_libs['dmumps'] = load_library('dmumps').dmumps_c
    mumps_libs['zmumps'] = load_library('zmumps').zmumps_c


def coo_is_symmetric(mtx, tol=1e-6):
    r, c = mtx.row, mtx.col
    odiag = nm.where(r != c)[0]
    if odiag.shape[0] == 0:
        return True

    odr, odc, odd = r[odiag], c[odiag], mtx.data[odiag]

    idxs = nm.where(odc > odr)[0]
    if (idxs.shape[0] * 2) == odiag.shape[0]:
        odr[idxs], odc[idxs] = odc[idxs], odr[idxs]

        idxs = nm.lexsort((odc, odr))
        odr, odc, odd = odr[idxs], odc[idxs], odd[idxs]
        d1 = nm.abs(nm.diff(nm.vstack([odr, odc]))).sum(axis=0)
        d2 = nm.diff(odd)
        if nm.all(d1[::2] == 0):
            vals = nm.abs(d2[::2])
            idxs = nm.where(vals > nm.finfo(vals.dtype).resolution)[0]
            odd_rs = odd.reshape((vals.shape[0], 2))
            vals[idxs] /= nm.abs(odd_rs).max(axis=1)[idxs]
            if nm.all(vals < tol):
                return True

    return False


class mumps_struc_c_x(ctypes.Structure):
    _fields_ = [
        ('sym', mumps_int),
        ('par', mumps_int),
        ('job', mumps_int),
        ('comm_fortran', mumps_int),
        ('icntl', mumps_int * 40),
        ('aux', ctypes.c_uint8 * AUX_LENGTH)
    ]


class mumps_struc_c_4(ctypes.Structure):  # MUMPS 4.10.0
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
        ('lsol_loc', mumps_int),
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
        ('version_number', ctypes.c_char * (14 + 1 + 1)),
        # /* For out-of-core */
        ('ooc_tmpdir', ctypes.c_char * 256),
        ('ooc_prefix', ctypes.c_char * 64),
        # /* To save the matrix in matrix market format */
        ('write_problem', ctypes.c_char * 256),
        ('lwk_user', mumps_int),
    ]


class mumps_struc_c_5_0(ctypes.Structure):  # MUMPS 5.0.x
    _fields_ = [
        ('sym', mumps_int),
        ('par', mumps_int),
        ('job', mumps_int),
        ('comm_fortran', mumps_int),
        ('icntl', mumps_int * 40),
        ('keep', mumps_int * 500),
        ('cntl', mumps_real * 15),
        ('dkeep', mumps_real * 130),
        ('keep8', mumps_int8 * 150),
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
        ('colsca_from_mumps', mumps_int),
        ('rowsca_from_mumps', mumps_int),
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
        ('lsol_loc', mumps_int),
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
        #  length=25 in FORTRAN + 1 for final \0 + 1 for alignment */
        ('version_number', ctypes.c_char * (25 + 1 + 1)),
        # /* For out-of-core */
        ('ooc_tmpdir', ctypes.c_char * 256),
        ('ooc_prefix', ctypes.c_char * 64),
        # /* To save the matrix in matrix market format */
        ('write_problem', ctypes.c_char * 256),
        ('lwk_user', mumps_int),
    ]


class mumps_struc_c_5_1(ctypes.Structure):  # MUMPS 5.1.x
    _fields_ = [
        ('sym', mumps_int),
        ('par', mumps_int),
        ('job', mumps_int),
        ('comm_fortran', mumps_int),
        ('icntl', mumps_int * 40),
        ('keep', mumps_int * 500),
        ('cntl', mumps_real * 15),
        ('dkeep', mumps_real * 230),
        ('keep8', mumps_int8 * 150),
        ('n', mumps_int),
        #
        ('nz_alloc', mumps_int),
        # /* Assembled entry */
        ('nz', mumps_int),
        ('nnz', mumps_int8),
        ('irn', mumps_pint),
        ('jcn', mumps_pint),
        ('a', mumps_pcomplex),
        # /* Distributed entry */
        ('nz_loc', mumps_int),
        ('nnz_loc', mumps_int8),
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
        ('colsca_from_mumps', mumps_int),
        ('rowsca_from_mumps', mumps_int),
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
        ('lsol_loc', mumps_int),
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
        #  length=30 in FORTRAN + 1 for final \0 + 1 for alignment */
        ('version_number', ctypes.c_char * (30 + 1 + 1)),
        # /* For out-of-core */
        ('ooc_tmpdir', ctypes.c_char * 256),
        ('ooc_prefix', ctypes.c_char * 64),
        # /* To save the matrix in matrix market format */
        ('write_problem', ctypes.c_char * 256),
        ('lwk_user', mumps_int),
        # /* For save/restore feature */
        ('save_dir', ctypes.c_char * 256),
        ('save_prefix', ctypes.c_char * 256)
    ]


class mumps_struc_c_5_2(ctypes.Structure):  # MUMPS 5.2.x
    _fields_ = [
        ('sym', mumps_int),
        ('par', mumps_int),
        ('job', mumps_int),
        ('comm_fortran', mumps_int),
        ('icntl', mumps_int * 60),
        ('keep', mumps_int * 500),
        ('cntl', mumps_real * 15),
        ('dkeep', mumps_real * 230),
        ('keep8', mumps_int8 * 150),
        ('n', mumps_int),
        #
        ('nz_alloc', mumps_int),
        # /* Assembled entry */
        ('nz', mumps_int),
        ('nnz', mumps_int8),
        ('irn', mumps_pint),
        ('jcn', mumps_pint),
        ('a', mumps_pcomplex),
        # /* Distributed entry */
        ('nz_loc', mumps_int),
        ('nnz_loc', mumps_int8),
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
        ('colsca_from_mumps', mumps_int),
        ('rowsca_from_mumps', mumps_int),
        # /* RHS, solution, ouptput data and statistics */
        ('rhs', mumps_pcomplex),
        ('redrhs', mumps_pcomplex),
        ('rhs_sparse', mumps_pcomplex),
        ('sol_loc', mumps_pcomplex),
        ('rhs_loc', mumps_pcomplex),
        ('irhs_sparse', mumps_pint),
        ('irhs_ptr', mumps_pint),
        ('isol_loc', mumps_pint),
        ('irhs_loc', mumps_pint),
        ('nrhs', mumps_int),
        ('lrhs', mumps_int),
        ('lredrhs', mumps_int),
        ('nz_rhs', mumps_int),
        ('lsol_loc', mumps_int),
        ('nloc_rhs', mumps_int),
        ('lrhs_loc', mumps_int),
        ('schur_mloc', mumps_int),
        ('schur_nloc', mumps_int),
        ('schur_lld', mumps_int),
        ('mblock', mumps_int),
        ('nblock', mumps_int),
        ('nprow', mumps_int),
        ('npcol', mumps_int),
        ('info', mumps_int * 80),
        ('infog', mumps_int * 80),
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
        #  length=30 in FORTRAN + 1 for final \0 + 1 for alignment */
        ('version_number', ctypes.c_char * (30 + 1 + 1)),
        # /* For out-of-core */
        ('ooc_tmpdir', ctypes.c_char * 256),
        ('ooc_prefix', ctypes.c_char * 64),
        # /* To save the matrix in matrix market format */
        ('write_problem', ctypes.c_char * 256),
        ('lwk_user', mumps_int),
        # /* For save/restore feature */
        ('save_dir', ctypes.c_char * 256),
        ('save_prefix', ctypes.c_char * 256),
        # /* Metis options */
        ('metis_options', mumps_int * 40),
    ]


class mumps_struc_c_5_3(ctypes.Structure):  # MUMPS 5.3.x
    _fields_ = [
        ('sym', mumps_int),
        ('par', mumps_int),
        ('job', mumps_int),
        ('comm_fortran', mumps_int),
        ('icntl', mumps_int * 60),
        ('keep', mumps_int * 500),
        ('cntl', mumps_real * 15),
        ('dkeep', mumps_real * 230),
        ('keep8', mumps_int8 * 150),
        ('n', mumps_int),
        ('nbl', mumps_int),
        #
        ('nz_alloc', mumps_int),
        # /* Assembled entry */
        ('nz', mumps_int),
        ('nnz', mumps_int8),
        ('irn', mumps_pint),
        ('jcn', mumps_pint),
        ('a', mumps_pcomplex),
        # /* Distributed entry */
        ('nz_loc', mumps_int),
        ('nnz_loc', mumps_int8),
        ('irn_loc', mumps_pint),
        ('jcn_loc', mumps_pint),
        ('a_loc', mumps_pcomplex),
        # /* Element entry */
        ('nelt', mumps_int),
        ('eltptr', mumps_pint),
        ('eltvar', mumps_pint),
        ('a_elt', mumps_pcomplex),
        # /* Matrix by blocks */
        ('blkptr', mumps_pint),
        ('blkvar', mumps_pint),
        # /* Ordering, if given by user */
        ('perm_in', mumps_pint),
        # /* Orderings returned to user */
        ('sym_perm', mumps_pint),
        ('uns_perm', mumps_pint),
        # /* Scaling (input only in this version) */
        ('colsca', mumps_preal),
        ('rowsca', mumps_preal),
        ('colsca_from_mumps', mumps_int),
        ('rowsca_from_mumps', mumps_int),
        # /* RHS, solution, ouptput data and statistics */
        ('rhs', mumps_pcomplex),
        ('redrhs', mumps_pcomplex),
        ('rhs_sparse', mumps_pcomplex),
        ('sol_loc', mumps_pcomplex),
        ('rhs_loc', mumps_pcomplex),
        ('irhs_sparse', mumps_pint),
        ('irhs_ptr', mumps_pint),
        ('isol_loc', mumps_pint),
        ('irhs_loc', mumps_pint),
        ('nrhs', mumps_int),
        ('lrhs', mumps_int),
        ('lredrhs', mumps_int),
        ('nz_rhs', mumps_int),
        ('lsol_loc', mumps_int),
        ('nloc_rhs', mumps_int),
        ('lrhs_loc', mumps_int),
        ('schur_mloc', mumps_int),
        ('schur_nloc', mumps_int),
        ('schur_lld', mumps_int),
        ('mblock', mumps_int),
        ('nblock', mumps_int),
        ('nprow', mumps_int),
        ('npcol', mumps_int),
        ('info', mumps_int * 80),
        ('infog', mumps_int * 80),
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
        #  length=30 in FORTRAN + 1 for final \0 + 1 for alignment */
        ('version_number', ctypes.c_char * (30 + 1 + 1)),
        # /* For out-of-core */
        ('ooc_tmpdir', ctypes.c_char * 256),
        ('ooc_prefix', ctypes.c_char * 64),
        # /* To save the matrix in matrix market format */
        ('write_problem', ctypes.c_char * 256),
        ('lwk_user', mumps_int),
        # /* For save/restore feature */
        ('save_dir', ctypes.c_char * 256),
        ('save_prefix', ctypes.c_char * 256),
        # /* Metis options */
        ('metis_options', mumps_int * 40),
    ]


class MumpsSolver(object):
    """MUMPS object."""

    def __init__(self, is_sym=False, mpi_comm=None,
                 system='real', silent=True, mem_relax=20):
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
        mem_relax : int
            The percentage increase in the estimated working space.
        """
        self.struct = None

        if not use_mpi:
            raise AttributeError('No mpi4py found! Required by MUMPS solver.')

        if len(mumps_libs) == 0:
            load_mumps_libraries()

        if system == 'real':
            self._mumps_c = mumps_libs['dmumps']
        elif system == 'complex':
            self._mumps_c = mumps_libs['zmumps']

        self.mpi_comm = MPI.COMM_WORLD if mpi_comm is None else mpi_comm
        self._mumps_c.restype = None

        # determine mumps version
        self._mumps_c.argtypes = [c_pointer(mumps_struc_c_x)]

        self.struct = mumps_struc_c_x()
        self.struct.par = 1
        self.struct.sym = 0
        self.struct.comm_fortran = self.mpi_comm.py2f()
        self.struct.job = -1

        self._mumps_c(ctypes.byref(self.struct))

        arr = nm.ctypeslib.as_array(self.struct.aux)
        idxs = nm.where(nm.logical_and(arr >= ord('.'), arr <= ord('9')))[0]
        s = dec(arr[idxs].tobytes())
        vnums = re.findall(r'^.*(\d)\.(\d+)\.\d+.*$', s)[-1]
        version = int(vnums[0]) + int(vnums[1]) * 1e-2
        if version < 5:
            mumps_struc_c = mumps_struc_c_4
        elif version >= 5 and version < 5.01:
            mumps_struc_c = mumps_struc_c_5_0
        elif version >= 5.01 and version < 5.02:
            mumps_struc_c = mumps_struc_c_5_1
        elif version >= 5.02 and version < 5.03:
            mumps_struc_c = mumps_struc_c_5_2
        elif version >= 5.03:
            mumps_struc_c = mumps_struc_c_5_3

        self.struct.job = -2

        self.set_silent()
        self._mumps_c(ctypes.byref(self.struct))

        # init mumps library
        self._mumps_c.argtypes = [c_pointer(mumps_struc_c)]

        self.struct = mumps_struc_c()
        self.struct.par = 1
        self.struct.sym = 2 if is_sym else 0
        self.struct.n = 0
        self.struct.comm_fortran = self.mpi_comm.py2f()

        self._mumps_call(job=-1)  # init

        self.rank = self.mpi_comm.rank
        self._data = {}

        # be silent
        if silent:
            self.set_silent()

        self.struct.icntl[13] = mem_relax

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

    def set_mtx_centralized(self, mtx):
        """
        Set the sparse matrix.

        Parameters
        ----------
        mtx : scipy sparse martix
            The sparse matrix in COO format.
        """
        assert mtx.shape[0] == mtx.shape[1]

        rr = mtx.row + 1
        cc = mtx.col + 1
        data = mtx.data

        if self.struct.sym > 0:
            idxs = nm.where(cc >= rr)[0]  # upper triangular matrix
            rr, cc, data = rr[idxs], cc[idxs], data[idxs]

        self.set_rcd_centralized(rr, cc, data, mtx.shape[0])

    def set_rcd_centralized(self, ir, ic, data, n):
        """
        Set the matrix by row and column indicies and data vector.
        The matrix shape is determined by the maximal values of
        row and column indicies. The indices start with 1.

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
        assert ir.shape[0] == ic.shape[0] == data.shape[0]

        self._data.update(ir=ir, ic=ic, data=data)
        self.struct.n = n
        self.struct.nz = ir.shape[0]
        if hasattr(self.struct, 'nnz'):
            self.struct.nnz = ir.shape[0]
        self.struct.irn = ir.ctypes.data_as(mumps_pint)
        self.struct.jcn = ic.ctypes.data_as(mumps_pint)
        self.struct.a = data.ctypes.data_as(mumps_pcomplex)

    def set_rhs(self, rhs):
        """Set the right hand side of the linear system."""
        self._data.update(rhs=rhs)
        self.struct.rhs = rhs.ctypes.data_as(mumps_pcomplex)
        self.struct.lrhs = rhs.shape[0]
        if len(rhs.shape) == 2:
            self.struct.nrhs = rhs.shape[1]

    def __call__(self, job):
        """Set the job and call MUMPS."""
        self._mumps_call(job)

    def get_schur(self, schur_list):
        """Get the Schur matrix and the condensed right-hand side vector.

        Parameters
        ----------
        schur_list : array
            The list of the Schur DOFs (indexing starts with 1).

        Returns
        -------
        schur_arr : array
            The Schur matrix of order 'schur_size'.
        schur_rhs : array
            The reduced right-hand side vector.
        """
        # Schur
        slist = schur_list + 1
        schur_size = slist.shape[0]
        schur_arr = nm.empty((schur_size**2, ), dtype='d')
        schur_rhs = nm.empty((schur_size, ), dtype='d')
        self._schur_rhs = schur_rhs

        self.struct.size_schur = schur_size
        self.struct.listvar_schur = slist.ctypes.data_as(mumps_pint)
        self.struct.schur = schur_arr.ctypes.data_as(mumps_pcomplex)
        self.struct.lredrhs = schur_size
        self.struct.redrhs = schur_rhs.ctypes.data_as(mumps_pcomplex)

        # get matrix
        self.struct.schur_lld = schur_size
        self.struct.nprow = 1
        self.struct.npcol = 1
        self.struct.mblock = 100
        self.struct.nblock = 100

        self.struct.icntl[18] = 3  # centr. Schur complement stored by columns
        self.struct.job = 4  # analyze + factorize
        self._mumps_c(ctypes.byref(self.struct))

        # get RHS
        self.struct.icntl[25] = 1  # Reduction/condensation phase
        self.struct.job = 3  # solve
        self._mumps_c(ctypes.byref(self.struct))

        return schur_arr.reshape((schur_size, schur_size)), schur_rhs

    def expand_schur(self, x2):
        """Expand the Schur local solution on the complete solution.

        Parameters
        ----------
        x2 : array
            The local Schur solution.

        Returns
        -------
        x : array
            The global solution.
        """
        self._schur_rhs[:] = x2
        self.struct.icntl[25] = 2  # Expansion phase
        self.struct.job = 3  # solve
        self._mumps_c(ctypes.byref(self.struct))

        return self._data['rhs']

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

        if ctypes is not None:
            self._mumps_c(ctypes.byref(self.struct))

        if self.struct.infog[0] < 0:
            raise RuntimeError("MUMPS error: %d" % self.struct.infog[0])

import numpy as nm

from sfepy.base.base import output, Struct
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.homogenization.homogen_app import HomogenizationApp
from sfepy.homogenization.coefficients import Coefficients
import tables as pt
from sfepy.discrete.fem.meshio import HDF5MeshIO
import sfepy.linalg as la
import sfepy.base.multiproc as multi
import os.path as op


def get_homog_coefs_linear(ts, coor, mode,
                           micro_filename=None, regenerate=False,
                           coefs_filename=None, define_args=None,
                           output_dir=None):

    oprefix = output.prefix
    output.prefix = 'micro:'

    required, other = get_standard_keywords()
    required.remove( 'equations' )

    conf = ProblemConf.from_file(micro_filename, required, other,
                                 verbose=False, define_args=define_args)
    if output_dir is not None:
        conf.options.output_dir = output_dir

    if coefs_filename is None:
        coefs_filename = conf.options.get('coefs_filename', 'coefs')
        coefs_filename = op.join(conf.options.get('output_dir', '.'),
                                 coefs_filename) + '.h5'

    if not regenerate:
        if op.exists( coefs_filename ):
            if not pt.is_hdf5_file( coefs_filename ):
                regenerate = True
        else:
            regenerate = True

    if regenerate:
        options = Struct( output_filename_trunk = None )

        app = HomogenizationApp( conf, options, 'micro:' )
        coefs = app()
        if type(coefs) is tuple:
            coefs = coefs[0]

        coefs.to_file_hdf5( coefs_filename )
    else:
        coefs = Coefficients.from_file_hdf5( coefs_filename )

    out = {}
    if mode == None:
        for key, val in coefs.__dict__.items():
            out[key] = val

    elif mode == 'qp':
        for key, val in coefs.__dict__.items():
            if type( val ) == nm.ndarray or type(val) == nm.float64:
                out[key] = nm.tile( val, (coor.shape[0], 1, 1) )
            elif type(val) == dict:
                for key2, val2 in val.items():
                    if type(val2) == nm.ndarray or type(val2) == nm.float64:
                        out[key+'_'+key2] = \
                                          nm.tile(val2, (coor.shape[0], 1, 1))

    else:
        out = None

    output.prefix = oprefix

    return out

def get_homog_coefs_nonlinear(ts, coor, mode, macro_data=None,
                              term=None, problem=None,
                              iteration=None, define_args=None,
                              output_dir=None, ret_corrs=False, **kwargs):
    if not (mode == 'qp'):
        return

    oprefix = output.prefix
    output.prefix = 'micro:'

    if not hasattr(problem, 'homogen_app'):
        required, other = get_standard_keywords()
        required.remove('equations')
        micro_file = problem.conf.options.micro_filename
        conf = ProblemConf.from_file(micro_file, required, other,
                                     verbose=False, define_args=define_args)
        if output_dir is not None:
            conf.options.output_dir = output_dir
        options = Struct(output_filename_trunk=None)
        app = HomogenizationApp(conf, options, 'micro:',
                                n_micro=coor.shape[0])
        problem.homogen_app = app

        if hasattr(app.app_options, 'use_mpi') and app.app_options.use_mpi:
            multiproc, multiproc_mode = multi.get_multiproc(mpi=True)
            multi_mpi = multiproc if multiproc_mode == 'mpi' else None
        else:
            multi_mpi = None

        app.multi_mpi = multi_mpi

        if multi_mpi is not None:
            multi_mpi.master_send_task('init', (micro_file, coor.shape[0]))
    else:
        app = problem.homogen_app
        multi_mpi = app.multi_mpi

    if macro_data is not None:
        macro_data['macro_time_step'] = ts.step

    app.setup_macro_data(macro_data)

    if multi_mpi is not None:
        multi_mpi.master_send_task('calculate', (macro_data, ts, iteration))

    coefs, deps = app(ret_all=True, itime=ts.step, iiter=iteration)

    if type(coefs) is tuple:
        coefs = coefs[0]

    out = {}
    for key, val in coefs.__dict__.items():
        if isinstance(val, list):
            out[key] = nm.array(val)
        elif isinstance(val, dict):
            for key2, val2 in val.items():
                out[key+'_'+key2] = nm.array(val2)

    for key in out.keys():
        shape = out[key].shape
        if len(shape) == 1:
            out[key] = out[key].reshape(shape + (1, 1))
        elif len(shape) == 2:
            out[key] = out[key].reshape(shape + (1,))

    output.prefix = oprefix

    if ret_corrs:
        return out, deps
    else:
        return out

def get_correctors_from_file_hdf5(coefs_filename='coefs.h5',
                                  dump_names=None):

    if dump_names == None:
        coefs = Coefficients.from_file_hdf5( coefs_filename )
        if hasattr(coefs, 'save_names'):
            dump_names = coefs.save_names
        else:
            raise ValueError( ' "filenames" coefficient must be used!' )

    out = {}

    for key, val in dump_names.items():
        if type(val) in [tuple, list]:
            h5name, corr_name = val
        else:
            h5name, corr_name = val, op.split(val)[-1]

        io = HDF5MeshIO(h5name + '.h5')
        try:
            ts = io.read_time_stepper()
        except ValueError:
            ts = None

        if ts is None:
            data = io.read_data(0)
            dkeys = list(data.keys())
            corr = {}
            for dk in dkeys:
                corr[dk] = data[dk].data.reshape(data[dk].shape)

            out[corr_name] = corr
        else:
            n_step = ts[3]
            out[corr_name] = []
            for step in range(n_step):
                data = io.read_data(step)
                dkeys = list(data.keys())
                corr = {}
                for dk in dkeys:
                    corr[dk] = data[dk].data.reshape(data[dk].shape)

                out[corr_name].append(corr)
            out[corr_name + '_ts'] = ts

    return out

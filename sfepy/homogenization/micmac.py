import numpy as nm

from sfepy.base.base import output, Struct
from sfepy.base.conf import ProblemConf, get_standard_keywords
from homogen import HomogenizationApp
from sfepy.homogenization.coefficients import Coefficients
import tables as pt
from sfepy.fem.meshio import HDF5MeshIO
import os.path as op

def get_homog_coefs_linear(ts, coor, mode,
                           micro_filename=None, regenerate=False,
                           coefs_filename=None):

    oprefix = output.prefix
    output.prefix = 'micro:'

    required, other = get_standard_keywords()
    required.remove( 'equations' )

    conf = ProblemConf.from_file(micro_filename, required, other, verbose=False)
    if coefs_filename is None:
        coefs_filename = conf.options.get_default_attr('coefs_filename', 'coefs')
        coefs_filename = op.join(conf.options.get_default_attr('output_dir', '.'),
                                 coefs_filename) + '.h5'

    if not regenerate:
        if op.exists( coefs_filename ):
            if not pt.isHDF5File( coefs_filename ):
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
        for key, val in coefs.__dict__.iteritems():
            out[key] = val

    elif mode == 'qp':
        for key, val in coefs.__dict__.iteritems():
            if type( val ) == nm.ndarray or type(val) == nm.float64:
                out[key] = nm.tile( val, (coor.shape[0], 1, 1) )
            elif type(val) == dict:
                for key2, val2 in val.iteritems():
                    if type(val2) == nm.ndarray or type(val2) == nm.float64:
                        out[key+'_'+key2] = \
                                          nm.tile(val2, (coor.shape[0], 1, 1))

    else:
        out = None

    output.prefix = oprefix

    return out

def get_correctors_from_file( coefs_filename = 'coefs.h5',
                              dump_names = None ):

    if dump_names == None:
        coefs = Coefficients.from_file_hdf5( coefs_filename )
        if hasattr( coefs, 'dump_names' ):
            dump_names = coefs.dump_names
        else:
            raise ValueError( ' "filenames" coefficient must be used!' )

    out = {}

    for key, val in dump_names.iteritems():
        corr_name = op.split( val )[-1]
        io = HDF5MeshIO( val+'.h5' )
        data = io.read_data( 0 )
        dkeys = data.keys()
        corr = {}
        for dk in dkeys:
            corr[dk] = data[dk].data.reshape(data[dk].shape)

        out[corr_name] = corr

    return out

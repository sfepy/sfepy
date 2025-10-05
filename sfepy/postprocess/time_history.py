import numpy as nm

from sfepy.base.base import output, OneTypeList, Struct
from sfepy.discrete.fem.mesh import Mesh
from sfepy.discrete.fem.meshio import MeshIO
from sfepy.solvers.ts import TimeStepper
from sfepy.base.ioutils import get_trunk, write_dict_hdf5

def _linearize(out, fields, linearization):
    new = {}
    for key, val in out.items():
        field = fields[val.field_name]
        new.update(field.create_output(val.data, var_name=key,
                                       dof_names=val.dofs, key=key,
                                       linearization=linearization))

    return new


def dump_to_vtk(filename, output_filename_trunk=None, step0=0, steps=None,
                fields=None, linearization=None):
    """Dump a multi-time-step results file into a sequence of VTK files."""
    def _save_step(suffix, out, mesh):
        if linearization is not None:
            output('linearizing...')
            out = _linearize(out, fields, linearization)
            output('...done')
            for key, val in out.items():
                lmesh = val.get('mesh', mesh)
                lmesh.write(output_filename_trunk + '_' + key + suffix,
                            io='auto', out={key : val})
                if hasattr(val, 'levels'):
                    output('max. refinement per group:', val.levels)

        else:
            mesh.write(output_filename_trunk + suffix, io='auto', out=out)

    output('dumping to VTK...')

    io = MeshIO.any_from_filename(filename)
    mesh = Mesh.from_file(filename, io=io)

    if output_filename_trunk is None:
        output_filename_trunk = get_trunk(filename)

    try:
        ts = TimeStepper(*io.read_time_stepper())
        all_steps, times, nts, dts = extract_times(filename)

    except ValueError:
        output('no time stepping info found, assuming single step')

        out = io.read_data(0)
        if out is not None:
            _save_step('.vtk', out, mesh)

        ret = None

    else:
        ts.times = times
        ts.n_step = times.shape[0]

        if steps is None:
            ii0 = nm.searchsorted(all_steps, step0)
            iterator = ((all_steps[ii], times[ii])
                        for ii in range(ii0, len(times)))

        else:
            iterator = [(step, ts.times[step]) for step in steps]

        max_step = all_steps.max()
        for step, time in iterator:
            output(ts.format % (step, max_step))
            out = io.read_data(step)
            if out is None: break

            _save_step('.' + ts.suffix % step + '.vtk', out, mesh)

        ret = ts.suffix

    output('...done')
    return ret

def extract_times(filename):
    """
    Read true time step data from individual time steps.

    Returns
    -------
    steps : array
        The time steps.
    times : array
        The times of the time steps.
    nts : array
        The normalized times of the time steps, in [0, 1].
    dts : array
        The true time deltas.
    """
    io = MeshIO.any_from_filename(filename)
    steps, times, nts = io.read_times()

    dts = nm.ediff1d(times, to_end=0)

    return steps, times, nts, dts

def extract_time_history(filename, extract, verbose=True):
    """Extract time history of a variable from a multi-time-step results file.

    Parameters
    ----------
    filename : str
        The name of file to extract from.
    extract : str
        The description of what to extract in a string of comma-separated
        description items. A description item consists of: name of the variable
        to extract, mode ('e' for elements, 'n' for nodes), ids of the nodes or
        elements (given by the mode). Example: 'u n 10 15, p e 0' means
        variable 'u' in nodes 10, 15 and variable 'p' in element 0.
    verbose : bool
        Verbosity control.

    Returns
    -------
    ths : dict
        The time histories in a dict with variable names as keys. If a nodal
        variable is requested in elements, its value is a dict of histories in
        the element nodes.
    ts : TimeStepper instance
        The time stepping information.
    """
    output('extracting selected data...', verbose=verbose)

    output('selection:', extract, verbose=verbose)

    ##
    # Parse extractions.
    pes = OneTypeList(Struct)
    for chunk in extract.split(','):
        aux = chunk.strip().split()
        pes.append(Struct(var=aux[0],
                          mode=aux[1],
                          indx=list(map(int, aux[2:]))))

    ##
    # Verify array limits.
    mesh = Mesh.from_file(filename)
    for pe in pes:
        if pe.mode == 'n':
            for ii in pe.indx:
                if (ii < 0) or (ii >= mesh.n_nod):
                    raise ValueError('node index 0 <= %d < %d!'
                                     % (ii, mesh.n_nod))

        if pe.mode == 'e':
            for ii, ie in enumerate(pe.indx[:]):
                if (ie < 0) or (ie >= mesh.n_el):
                    raise ValueError('element index 0 <= %d < %d!'
                                     % (ie, mesh.n_el))
                pe.indx[ii] = ie

    ##
    # Extract data.
    io = MeshIO.any_from_filename(filename)
    ths = {}
    for pe in pes:
        mode, nname = io.read_data_header(pe.var)
        output(mode, nname, verbose=verbose)

        if ((pe.mode == 'n' and mode == 'vertex') or
            (pe.mode == 'e' and mode == 'cell')):
            th = io.read_time_history(nname, pe.indx)

        elif pe.mode == 'e' and mode == 'vertex':
            conn = mesh.conns[0]
            th = {}
            for iel in pe.indx:
                ips = conn[iel]
                th[iel] = io.read_time_history(nname, ips)
        else:
            raise ValueError('cannot extract cell data %s in nodes!' % pe.var)

        ths[pe.var] = th

    output('...done', verbose=verbose)

    ts = TimeStepper(*io.read_time_stepper())
    # Force actual times.
    steps, times, nts, dts = extract_times(filename)
    ts.times = times
    ts.nt = nts

    return ths, ts

def average_vertex_var_in_cells(ths_in):
    """Average histories in the element nodes for each nodal variable
        originally requested in elements."""
    ths = dict.fromkeys(list(ths_in.keys()))
    for var, th in ths_in.items():
        aux = dict.fromkeys(list(th.keys()))
        for ir, data in th.items():
            if isinstance(data, dict):
                for ic, ndata in data.items():
                    if aux[ir] is None:
                        aux[ir] = ndata
                    else:
                        aux[ir] += ndata
                aux[ir] /= float(len(data))
            else:
                aux[ir] = data
        ths[var] = aux

    return ths

def save_time_history(ths, ts, filename_out):
    """Save time history and time-stepping information in a HDF5 file."""
    ths.update({'times' : ts.times, 'dt' : ts.dt})
    write_dict_hdf5(filename_out, ths)

def guess_time_units(times):
    """
    Given a vector of times in seconds, return suitable time units and
    new vector of times suitable for plotting.

    Parameters
    ----------
    times : array
        The vector of times in seconds.

    Returns
    -------
    new_times : array
        The vector of times in `units`.
    units : str
        The time units.
    """
    times = nm.asarray(times)

    if (times[-1] / 60.0 / 60.0) > 10.0:
        units = 'hours'
        new_times = times / 60.0 / 60.0

    elif (times[-1] / 60.0) > 10.0:
        units = 'min.'
        new_times = times / 60.0

    else:
        units = 's'
        new_times = times

    return new_times, units

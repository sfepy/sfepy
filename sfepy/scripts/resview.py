#!/usr/bin/env python
"""
This is a script for quick VTK-based visualizations of finite element
computations results.

In the examples below it is supposed that sfepy is installed. When using the
in-place build, replace ``sfepy-view`` by ``python3 sfepy/scripts/resview.py``.

Examples
--------
The examples assume that
``python -c "import sfepy; sfepy.test('--output-dir=output-tests')"``
has been run successfully and the resulting data files are present.

- View data in output-tests/test_navier_stokes.vtk::

    sfepy-view output-tests/navier_stokes-navier_stokes.vtk

- Customize the above output:
  plot0: field "p", switch on edges,
  plot1: field "u", surface with opacity 0.4, glyphs scaled by factor 2e-2::

    sfepy-view output-tests/navier_stokes-navier_stokes.vtk -f p:e:p0 u:o.4:p1 u:g:f2e-2:p1

- As above, but glyphs are scaled by the factor determined automatically as
  20% of the minimum bounding box size::

    sfepy-view output-tests/navier_stokes-navier_stokes.vtk -f p:e:p0 u:o.4:p1 u:g:f10%:p1

- View data and take a screenshot::

    sfepy-view output-tests/diffusion-poisson.vtk -o image.png

- Take a screenshot without a window popping up::

    sfepy-view output-tests/diffusion-poisson.vtk -o image.png --off-screen

- Create animation from output-tests/diffusion-time_poisson.*.vtk::

    sfepy-view output-tests/diffusion-time_poisson.*.vtk -a mov.mp4

- Create animation from output-tests/test_hyperelastic.*.vtk,
  set frame rate to 3, plot displacements and mooney_rivlin_stress::

    sfepy-view output-tests/test_hyperelastic_TL.*.vtk -f u:wu:e:p0 mooney_rivlin_stress:p1 -a mov.mp4 -r 3
"""
from argparse import ArgumentParser, Action, RawDescriptionHelpFormatter
from ast import literal_eval
import numpy as nm
import os.path as osp

import pyvista as pv
from vtk.util.numpy_support import numpy_to_vtk

cache = {}


def get_camera_position(bounds, azimuth, elevation, distance=None, zoom=1.):
    phi, psi = nm.deg2rad(azimuth), nm.deg2rad(elevation)
    bounds = nm.asarray(bounds)

    if distance is not None:
        r = distance / zoom
    else:
        r = max(bounds[1::2] - bounds[::2]) * 2.0 / zoom

    center = (bounds[1::2] + bounds[::2]) * 0.5

    # camera position
    position = (r * nm.cos(phi) * nm.sin(psi),
                r * nm.sin(phi) * nm.sin(psi),
                r * nm.cos(psi))

    # view up
    view_up = (0, 0, 1)
    if abs(elevation) < 5. or abs(elevation) > 175.:
        view_up = (nm.sin(phi), nm.cos(phi), 0)

    return [position, tuple(center), view_up]


def parse_options(opts, separator=':'):
    out = {}
    if opts is None:
        return out

    for v in opts.split(separator):
        if len(v) < 2:
            val = True
        elif v[1:].isalpha():
            val = v[1:]
        elif v[-1] == '%':
            val = ('%', float(v[1:-1]))
        else:
            try:
                val = literal_eval(v[1:])
            except ValueError:
                val = v[1:]

        out[v[0]] = val

    return out


def make_title(filenames):
    title = ', '.join(filenames)
    title = title if len(title) < 80 else title[:77] + '...'
    return title


def make_cells_from_conn(conns, convert_to_vtk_type):
    cells, cell_type, offset = [], [], []
    _offset = 0
    for ctype, conn in conns.items():
        nc, np = conn.shape

        aux = nm.empty((nc, np + 1), dtype=int)
        aux[:, 0] = np
        aux[:, 1:] = conn
        cells.append(aux.ravel())

        cell_type.append(nm.full(nc, convert_to_vtk_type[ctype]))
        offset.append(nm.arange(nc) * (np + 1) + _offset)
        _offset += nc

    cells = nm.concatenate(cells)
    cell_type = nm.concatenate(cell_type)
    offset = nm.concatenate(offset)

    return cells, cell_type, offset


def add_mat_id_to_grid(grid, cell_groups):
    val = numpy_to_vtk(cell_groups)
    val.SetName('mat_id')
    grid.GetCellData().AddArray(val)
    return grid


vtk_cell_types = {'1_1': 1, '1_2': 3, '2_2': 3, '3_2': 3,
                  '2_3': 5, '2_4': 9, '3_4': 10, '3_8': 12}

def make_grid_from_mesh(mesh, add_mat_id=False):
    desc = mesh.descs[0]
    nv, dim = mesh.coors.shape

    points = nm.c_[mesh.coors, nm.zeros((nv, 3 - dim))]
    cells, cell_type, offset = make_cells_from_conn(
        {desc: mesh.get_conn(desc)}, vtk_cell_types,
    )
    try:
        grid = pv.UnstructuredGrid(offset, cells, cell_type, points)

    except TypeError: # Pyvista >= 0.39.0
        grid = pv.UnstructuredGrid(cells, cell_type, points)

    if add_mat_id:
        add_mat_id_to_grid(grid, mesh.cmesh.cell_groups)

    return grid

def read_mesh(filenames, step=None, print_info=True, ret_n_steps=False,
              use_cache=True):
    _, ext = osp.splitext(filenames[0])
    if ext in ['.vtk', '.vtu']:
        fstep = 0 if step is None else step
        fname = filenames[fstep]
        key = (fname, fstep)
        if key not in cache or not use_cache:
            cache[key] = ((fstep, float(fstep)), pv.UnstructuredGrid(fname))
        ftime, mesh = cache[key]
        cache['n_steps'] = len(filenames)
    elif ext in ['.xdmf', '.xdmf3']:
        import meshio
        try:
            from meshio._common import meshio_to_vtk_type

        except ImportError:
            from meshio._vtk_common import meshio_to_vtk_type

        fname = filenames[0]
        key = (fname, step)
        if key not in cache:
            reader = meshio.xdmf.TimeSeriesReader(fname)
            points, _cells = reader.read_points_cells()
            points = nm.asarray(points)
            if points.shape[1] < 3:
                points = nm.pad(points, [(0, 0), (0, 3 - points.shape[1])])
            _dcells = {ct.type: ct.data for ct in _cells}

            cells, cell_type, offset = make_cells_from_conn(
                _dcells, meshio_to_vtk_type,
            )

            if not reader.num_steps:
                grid = pv.UnstructuredGrid(offset, cells, cell_type, points)
                add_mat_id_to_grid(grid, mesh.cmesh.cell_groups)
                cache[(fname, 0)] = (0.0, grid)

            grids = {}
            time = []
            for _step in range(reader.num_steps):
                grid = pv.UnstructuredGrid(offset, cells, cell_type, points)
                t, pd, cd = reader.read_data(_step)
                for dk, dv in pd.items():
                    val = numpy_to_vtk(dv)
                    val.SetName(dk)
                    grid.GetPointData().AddArray(val)

                for dk, dv in cd.items():
                    val = numpy_to_vtk(nm.vstack(dv).squeeze())
                    val.SetName(dk)
                    grid.GetCellData().AddArray(val)

                grids[t] = grid
                time.append(t)

            time.sort()
            for _step, t in enumerate(time):
                cache[(fname, _step)] = ((_step, t), grids[t])

            cache[(fname, None)] = cache[(fname, 0)]
            cache['n_steps'] = reader.num_steps

        ftime, mesh = cache[key]

    elif ext in ['.h5', '.h5x']:
        # Custom sfepy format.
        fname = filenames[0]
        if 'io' not in cache:
            from sfepy.discrete.fem.meshio import MeshIO

            io = MeshIO.any_from_filename(fname)

            smesh = io.read()

            cache['steps'], cache['times'], _ = io.read_times()
            if not len(cache['steps']):
                grid0 = make_grid_from_mesh(smesh, add_mat_id=True)
                cache[(fname, 0)] = ((0, 0.0), grid0)

            else:
                grid0 = make_grid_from_mesh(smesh, add_mat_id=False)

            cache['io'] = io
            cache['grid0'] = grid0
            cache['n_steps'] = len(cache['steps'])

        if step is None:
            step = 0
        key = (fname, step)
        if key not in cache:
            io = cache['io']
            grid = cache['grid0'].copy()

            _step = cache['steps'][step]
            datas = io.read_data(_step)
            for dk, data in datas.items():
                vval = data.data
                if 1 < len(data.dofs) < 3:
                    vval = nm.c_[vval,
                                 nm.zeros((len(vval), 3 - len(data.dofs)))]

                if data.mode == 'vertex':
                    val = numpy_to_vtk(vval)
                    val.SetName(dk)
                    grid.GetPointData().AddArray(val)

                else:
                    val = numpy_to_vtk(vval[:, 0, :, 0])
                    val.SetName(dk)
                    grid.GetCellData().AddArray(val)

            cache[(fname, step)] = ((_step, cache['times'][step]), grid)

        ftime, mesh = cache[key]

    else:
        fname = filenames[0]
        key = (fname, step)
        if key not in cache:
            from sfepy.discrete.fem.meshio import MeshIO
            from sfepy.discrete.fem import Mesh

            io = MeshIO.any_from_filename(fname)
            smesh = Mesh(fname)
            smesh = io.read(smesh)

            grid = make_grid_from_mesh(smesh, add_mat_id=True)
            cache[(fname, 0)] = ((0, 0.0), grid)
            cache['n_steps'] = len(filenames)

        ftime, mesh = cache[key]

    if print_info:
        arrs = {'s': [], 'v': [], 'o': []}
        for aname in mesh.array_names:
            if len(mesh[aname].shape) == 1 or mesh[aname].shape[1] == 1:
                arrs['s'].append(aname)
            elif mesh[aname].shape[1] == 3:
                arrs['v'].append(aname)
            else:
                arrs['o'].append(aname + '(%d)' % mesh[aname].shape[1])

        step_info = ' (step %d)' % step if step else ''
        print('mesh from %s%s:' % (fname, step_info))
        print('  points:  %d' % mesh.n_points)
        print('  cells:   %d' % mesh.n_cells)
        print('  bounds:  %s' % list(zip(nm.min(mesh.points, axis=0),
                                         nm.max(mesh.points, axis=0))))
        if len(arrs['s']) > 0:
            print('  scalars: %s' % ', '.join(arrs['s']))
        if len(arrs['v']) > 0:
            print('  vectors: %s' % ', '.join(arrs['v']))
        if len(arrs['o']) > 0:
            print('  others:  %s' % ', '.join(arrs['o']))
        print('  steps:   %d' % cache['n_steps'])

    if ret_n_steps:
        return ftime, mesh, cache['n_steps']
    else:
        return ftime, mesh

def ensure3d(arr):
    arr = nm.asanyarray(arr)
    return nm.pad(arr, (arr.ndim - 1) * [(0, 0)] + [(0, 3 - arr.shape[-1])])

def get_nonzero_norm(arr, **kwargs):
    norm = nm.linalg.norm(arr, **kwargs)
    if isinstance(norm, nm.ndarray):
        norm = norm.max()

    return norm if norm != 0 else 1.0

def read_probes_as_annotations(filenames, add_label=True):
    from sfepy.discrete.probes import read_results
    from sfepy.linalg.geometry import get_perpendiculars

    annotations = []
    for filename in filenames:
        header, results = read_results(filename)
        data = header.details
        if header.probe_class == 'PointsProbe':
            data = ensure3d(data)
            ann = [('points', data)]
            tcoors = data.mean(axis=0)

        elif header.probe_class == 'LineProbe':
            data = ensure3d(data)
            ann = [('line', data[0], data[1])]
            tcoors = data.mean(axis=0)

        elif header.probe_class == 'RayProbe':
            data[:2] = ensure3d(data[:2])
            ann = [('arrow', data[0], data[1])]
            if data[2]:
                ann += [('arrow', data[0], -data[1])]
            tcoors = data[0]

        elif header.probe_class == 'CircleProbe':
            data[:2] = ensure3d(data[:2])
            ann = [('disc', data[0], data[1], data[2])]
            vec = get_perpendiculars(data[1])[0]
            tcoors = data[0] + data[2] * vec

        else:
            raise ValueError(f'unknown probe kind! {header.probe_class}')

        if add_label:
            ann += [('text', [tcoors],
                     [osp.splitext(osp.basename(filename))[0]])]

        annotations.extend(ann)

    return annotations

def make_glyphs(mesh, scalars_name=None, geom_name=None):
    if (geom_name is None) or (geom_name == 'arrows'):
        glyphs = mesh.arrows

    else:
        vectors, vectors_name = mesh.active_vectors, mesh.active_vectors_name
        if vectors is None or vectors_name is None:
            return None

        if vectors.ndim != 2:
            raise ValueError('active vectors are not vectors.')

        if scalars_name is None:
            scalars_name = f'{vectors_name} Magnitude'
            scale = nm.linalg.norm(vectors, axis=1)
            mesh.point_data.set_array(scale, scalars_name)

        if geom_name == 'cylinders':
            geom = pv.Cylinder(
                direction=(1, 0, 0), height=1.0, radius=0.05, resolution=8,
            )

        elif geom_name == 'cones':
            geom = pv.Cone(
                direction=(1, 0, 0), height=1.0, radius=0.05, resolution=8,
            )
        elif geom_name == 'lines':
            geom = pv.Line(pointa=(0, 0, -0.5), pointb=(0, 0, 0.5))

        elif geom_name == 'tubes':
            line = pv.Line(pointa=(-0.5, 0, 0), pointb=(0.5, 0, 0))
            geom = line.tube(radius=0.05, n_sides=8)

        elif geom_name == 'spheres':
            geom = pv.Sphere(radius=0.5, theta_resolution=8, phi_resolution=8)

        else:
            raise ValueError(f'unsupported glyph geometry! {geom_name}')

        glyphs = mesh.glyph(orient=vectors_name, scale=scalars_name, geom=geom)

    return glyphs

def pv_plot(filenames, options, plotter=None, step=None, annotations=None,
            scalar_bar_limits=None, ret_scalar_bar_limits=False,
            step_inc=None, use_cache=True):
    fstep = (step if step is not None else options.step)
    if step_inc is not None:
        fstep += step_inc
    if fstep < 0:
        return
    if hasattr(plotter, 'resview_n_steps'): # Works for None as well.
        if fstep >= plotter.resview_n_steps:
            return

    plots = {}
    color = None

    if plotter is None:
        plotter = pv.Plotter(title=make_title(filenames))

    if step_inc is not None:
        plotter.clear()

    ftime, mesh, n_steps = read_mesh(filenames, fstep, ret_n_steps=True,
                                     use_cache=use_cache)
    steps = {fstep: mesh}
    ftimes = {fstep: ftime}

    bbox_sizes = nm.diff(nm.reshape(mesh.bounds, (-1, 2)), axis=1)
    dim = len(bbox_sizes)
    ii = nm.where(bbox_sizes > 0)[0]
    tdim = len(ii)
    if tdim == 0:
        ipv2, ipv = 1, 2
        print('WARNING: zero size mesh!')

    elif tdim > 1:
        ipv2, ipv = ii[-2:]

    else:
        ipv2, ipv = 0, 1

    if options.grid_vector1 is None:
        options.grid_vector1 = [0, 0, 0]
        options.grid_vector1[ipv] = 1.6

    if options.grid_vector2 is None:
        options.grid_vector2 = [0, 0, 0]
        options.grid_vector2[ipv2] = 1.6

    plotter.resview_step, plotter.resview_n_steps = fstep, n_steps

    fields_map = {}
    if len(options.fields_map) > 0:
        for cg, fields in options.fields_map:
            for field in fields.split(','):
                fields_map[field.strip()] = int(cg)

    if len(options.fields) == 0:
        fields = []
        position = 0
        for field in steps[fstep].array_names:
            if field in ['node_groups', 'mat_id']:
                continue

            fval = steps[fstep][field]
            bnds = steps[fstep].bounds
            mesh_size = (nm.array(bnds[1::2]) - nm.array(bnds[::2])).max()
            is_vector_field = (len(fval.shape) == 2) and (fval.shape[1] == dim)
            is_point_field = fval.shape[0] == steps[fstep].n_points
            if is_vector_field and is_point_field:
                scale = mesh_size * 0.15 / get_nonzero_norm(fval, axis=1)
                if not nm.isfinite(scale):
                    scale = 1.0
                fields.append((field, 'vs:o.4:p%d' % position))
                fields.append((field, 'g:f%e:p%d' % (scale, position)))
            else:
                fields.append((field, 'p%d' % position))

            position += 1

        if (len(fields) == 0):
            if  'mat_id' in mesh.array_names:
                fields.append(('mat_id', 'p0'))

            else:
                fields.append(('0', 'p0'))

    else:
        fields = options.fields

    plot_id = 0

    scalar_bars = {}
    for field, fopts in fields:
        opts = parse_options(fopts)
        plot_info = []

        if field == '0':
            field = None
            color = 'white'

        if field == '1':
            field = None
            color = 'black'

        if 's' in opts and step is None:  # plot data from a given step
            fstep = opts['s']

        if fstep not in steps:
            ftimes[fstep], steps[fstep] = read_mesh(filenames, step=fstep,
                                                    use_cache=use_cache)

        pipe = [steps[fstep].copy()]
        pos_bnds = pipe[0].bounds

        if field in fields_map:  # subregion
            mat_val = fields_map[field]
        elif 'm' in opts:
            mat_val = opts['m']
        else:
            mat_val = None

        if mat_val:
            if isinstance(mat_val, int):
                mat_val = [mat_val, mat_val]

            pipe.append(pipe[-1].threshold(value=mat_val,
                        scalars='mat_id', preference='cell'))

        if 'r' in opts:  # recalculate cell data to point data
            pipe.append(pipe[-1].cell_data_to_point_data())

        opacity = opts.get('o', options.opacity)  # mesh opacity
        show_edges = opts.get('e', options.show_edges)  # edge visibility
        style = {'s': 'surface',
                 'w': 'wireframe',
                 'p': 'points'}[opts.get('v', 's')]  # set style

        warp = opts.get('w', options.warp)  # warp mesh
        factor = opts.get('f', options.factor)
        if isinstance(factor, tuple):
            ws = nm.diff(nm.reshape(pipe[-1].bounds, (-1, 2)), axis=1)
            size = ws[ws > 0.0].min()
            factor_field = warp if warp is not None else field
            fmax = get_nonzero_norm(pipe[-1][factor_field], ord=nm.inf)
            factor = 0.01 * float(factor[1]) * size / fmax

        if warp:
            field_data = pipe[-1][warp]
            if field_data.ndim == 1:
                field_data.shape = (-1, 1)
            nc = field_data.shape[1]
            if nc == 1:  # Warp by scalar.
                pipe.append(pipe[-1].copy())
                pipe[-1].points[:, 2] += field_data[:, 0] * factor

            elif nc == 3:
                pipe.append(pipe[-1].copy())
                pipe[-1].points += field_data * factor

            else:
                raise ValueError('warp mesh: scalar or vector field required!')

            plot_info.append('warp=%s, factor=%.2e' % (warp, factor))

        position = opts.get('p', 0)  # determine plotting slot
        if 'p' in opts:
            size = nm.array(pos_bnds[1::2]) - nm.array(pos_bnds[::2])
            size = nm.maximum(size, 1e-2 * size.max())
            pipe.append(pipe[-1].copy())
            pos1 = position % options.max_plots
            pos2 = position // options.max_plots
            shift = pos1 * size * nm.array(options.grid_vector1)
            shift += pos2 * size * nm.array(options.grid_vector2)
            pipe[-1].translate(shift, inplace=True)

        if opts.get('l', options.outline):  # outline
            plotter.add_mesh(pipe[-1].outline(), color='k')

        is_vector_field = ((field is not None)
                           and (len(pipe[-1][field].shape) > 1)
                           and (pipe[-1][field].shape[1] == dim))
        is_point_field = (field is not None and
                          pipe[-1][field].shape[0] == pipe[-1].n_points)

        has_components = ((field is not None)
                          and (len(pipe[-1][field].shape) > 1)
                          and (pipe[-1][field].shape[1] > 1))
        if has_components:
            field_data = pipe[-1][field]
            scalar = field + '_magnitude'
            scalar_label = f'{field}:mag'
            pipe[-1][scalar] = nm.linalg.norm(field_data, axis=1)

        else:
            scalar = field
            scalar_label = scalar

        if 'g' in opts and is_point_field:  # glyphs
            geom_name, scalar_name = None, None
            if isinstance(opts['g'], str):
                aux = opts['g'].split(',')
                gfield = aux[0]
                if len(aux) > 2:
                    # Scaling scalar name.
                    scalar_name = aux[2]
                if len(aux) > 1:
                    # Glyph geometry name name.
                    geom_name = aux[1]

            else:
                gfield = opts['g']

            if isinstance(gfield, str):
                is_gvector_field = ((gfield is not None)
                                    and (len(pipe[-1][gfield].shape) > 1)
                                    and (pipe[-1][gfield].shape[1] == dim))

            else:
                gfield = field
                is_gvector_field = is_vector_field

            if is_gvector_field:
                pipe[-1][gfield] *= factor
                pipe[-1].set_active_vectors(gfield)

                pipe.append(make_glyphs(
                    pipe[-1],
                    scalars_name=scalar_name,
                    geom_name=geom_name,
                ))
                show_edges = False
                plot_info.append('glyphs=%s, factor=%.2e' % (field, factor))

            else:
                g_pipe = pipe[-1].compute_derivative(scalars=gfield)
                g_pipe['gradient'] *= factor
                g_pipe.set_active_vectors('gradient')
                pipe.append(make_glyphs(
                    g_pipe,
                    scalars_name=scalar_name,
                    geom_name=geom_name,
                ))
                show_edges = False
                plot_info.append('glyphs=grad(%s), factor=%.2e'
                                 % (field, factor))

                scalar = 'gradient Magnitude'
                scalar_label = f'grad({field}):mag'

        elif 'c' in opts and has_components:  # select field component
            comp = opts['c']
            scalar = field + '_%d' % comp
            pipe[-1][scalar] = field_data[:, comp]
        elif 't' in opts:  # streamlines
            scalar = field
            tube_radius = None
            opt = opts['t']
            if opt is True:
                npts = 20

            elif isinstance(opt, int):
                npts = opt

            else:
                aux = opt.split(',')
                npts = int(aux[0])
                if len(aux) > 2:
                    tube_radius = float(aux[2])
                if len(aux) > 1:
                    scalar = aux[1]

            if is_vector_field:
                sl_vector = field
                sl_pipe = pipe[-1]
            else:
                sl_vector = 'gradient'
                sl_pipe = pipe[-1].compute_derivative(scalars=field)

            cmin, cmax = sl_pipe.bounds[::2], sl_pipe.bounds[1::2]
            if tdim == 2:
                streamlines = sl_pipe.streamlines(vectors=sl_vector,
                                                  pointa=cmin, pointb=cmax,
                                                  n_points=npts,
                                                  max_time=1e12)

            else:
                radius = 0.5 * nm.linalg.norm(nm.array(cmax) - nm.array(cmin))
                streamlines = sl_pipe.streamlines(vectors=sl_vector,
                                                  source_radius=radius,
                                                  n_points=npts,
                                                  max_time=1e12)

            if tube_radius is None:
                pipe.append(streamlines)

            else:
                pipe.append(streamlines.tube(radius=tube_radius))

        isosurfaces = int(opts.get('i', options.isosurfaces))
        if isosurfaces > 0:  # iso-surfaces
            pipe[-1].set_active_scalars(scalar)
            field_data = pipe[-1][scalar]
            pars = (nm.min(field_data), nm.max(field_data), isosurfaces + 1)
            pipe.append(pipe[-1].contour(nm.linspace(*pars)))

        kwargs = {}
        if options.color_limits is not None:
            kwargs['clim'] = options.color_limits

        plotter.add_mesh(pipe[-1], scalars=scalar, color=color,
                         style=style, show_edges=show_edges,
                         opacity=opacity,
                         cmap=options.color_map,
                         show_scalar_bar=False, label=scalar_label,
                         **kwargs)

        bnds = pipe[-1].bounds
        if position not in plots:
            plots[position] = []

        plot_info = ':' + ','.join(plot_info) if len(plot_info) > 0 else ''
        plot_info = '%s(step %d)%s' % (scalar, fstep, plot_info)
        plots[position].append(((bnds[::2], bnds[1::2]), plot_info))

        if options.show_scalar_bars and scalar:
            if scalar not in scalar_bars:
                scalar_bars[scalar_label] = []

            field_data = pipe[-1][scalar]
            limits = (nm.min(field_data), nm.max(field_data))
            scalar_bars[scalar_label].append((limits, plotter.mapper,
                                              position))

        plot_id += 1

    if options.show_scalar_bars:
        if scalar_bar_limits is None:
            scalar_bar_limits = {}

            if options.color_limits is not None:
                scalar_bar_limits = {k: options.color_limits
                                     for k in scalar_bars.keys()}

        width, height = options.scalar_bar_size
        position_x, position_y, shift_x, shift_y = options.scalar_bar_position
        nslots = len(scalar_bars)
        for k, vs in scalar_bars.items():
            clim = scalar_bar_limits.get(k, (None, None))
            if clim[0] is None:
               clim = (nm.min([v[0] for v, _, _ in vs]), clim[1])
            if clim[1] is None:
               clim = (clim[0], nm.max([v[1] for v, _, _ in vs]))

            for _, mapper, _ in vs:
                mapper.scalar_range = clim
            _, mapper, slot = vs[0]

            slot_x = (nslots - slot - 1) if shift_x < 0 else slot
            x_pos = position_x + slot_x * width * shift_x
            slot_y = (nslots - slot - 1) if shift_y < 0 else slot
            y_pos = position_y + slot_y * height * shift_y

            plotter.add_scalar_bar(title=k,
                                   position_x=x_pos, position_y=y_pos,
                                   width=width, height=height,
                                   n_labels=2, mapper=mapper)
    if annotations is not None:
        for ann in annotations:
            kind, data = ann[0], ann[1:]
            if kind == 'points':
                pdata = pv.PolyData(data[0])

            elif kind == 'line':
                pdata = pv.Line(*data)

            elif kind == 'arrow':
                # Scale arrows by bbox.
                vec = data[1]
                maxs = bbox_sizes.max()
                vec *= 0.1 * maxs / get_nonzero_norm(vec)
                pdata = pv.Arrow(data[0], vec, scale='auto')

            elif kind == 'disc':
                pdata = pv.Disc(data[0],
                                inner=0.99 * data[2],
                                outer=1.01 * data[2],
                                normal=data[1], c_res=20)
            elif kind == 'text':
                plotter.add_point_labels(data[0], data[1],
                                         always_visible=True)

            else:
                raise ValueError(f'unknown annotation! {kind}')

            plotter.add_mesh(pdata, opacity=opacity, color='black',
                             line_width=5, render_lines_as_tubes=True,
                             point_size=10, render_points_as_spheres=True,
                             show_scalar_bar=False)

    if options.show_labels and len(plots) > 1:
        labels, points = [], []
        for k, v in plots.items():
            bnds = (nm.min(nm.array([iv[0][0] for iv in v]), axis=0),
                    nm.max(nm.array([iv[0][1] for iv in v]), axis=0))
            labels.append('plot:%d' % k)
            size = bnds[1] - bnds[0]
            olpos = options.label_position
            points.append(bnds[0] + nm.array(olpos[:3]) * size * olpos[3])

        plotter.add_point_labels(nm.array(points), labels)

    if options.show_step_time:
        _step, _time = ftimes[fstep]
        step_info = f'{_step:6d}: {_time:.8e}'
        if len(filenames) > 1:
            step_info += ': ' + osp.splitext(osp.basename(filenames[fstep]))[0]

        plotter.add_text(step_info, font='courier', font_size=10)

    for k, v in plots.items():
        print('plot %d: %s' % (k, '; '.join(iv[1] for iv in v)))

    if ret_scalar_bar_limits:
        return plotter, scalar_bar_limits
    else:
        return plotter


def print_camera_position(plotter):
    cp = nm.array([k for k in plotter.camera_position]).ravel()
    cp = ','.join(['%g' % k for k in cp])
    print(f'--camera-position="{cp}"')


def _get_cpos(plotter, options, camera_default=(225, 75, 0.9)):
    """
    Uses `plotter.bounds`, so call only after adding all meshes to the plotter.
    """
    if options.camera_position is not None:
        cpos = nm.array(options.camera_position)
        cpos = cpos.reshape((3, 3))
    elif options.camera:
        zoom = options.camera[2] if len(options.camera) > 2 else 1.
        cpos = get_camera_position(plotter.bounds,
                                   options.camera[0],
                                   options.camera[1],
                                   zoom=zoom)
    elif options.view_2d:
        cpos = None
    else:
        cpos = get_camera_position(plotter.bounds, camera_default[0],
                                   camera_default[1], zoom=camera_default[2])

    return cpos


class OptsToListAction(Action):
    separator = '='

    def __call__(self, parser, namespace, values, option_string=None):
        out = []
        for item in values:
            s = item.split(self.separator, 1)
            out.append((s[0].strip(), s[1].strip() if len(s) > 1 else None))

        setattr(namespace, self.dest, out)


class FieldOptsToListAction(OptsToListAction):
    separator = ':'


class StoreNumberAction(Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, literal_eval(values))


helps = {
    'fields':
        """fields to plot, options separated by ':' are possible:\n
        'g[,geom[,X]]' - plot glyphs (arrows by default) for vector or scalar
                         fields, where geom is one of ['arrows', 'cylinders',
                         'cones', 'lines', 'tubes', 'spheres'], scale by factor
                         F, and optionally by X scalar;
        'iX' - plot X isosurfaces;
        'tX[,Y]' - plot X streamlines, color optionally by Y, gradient employed
                   for scalar fields X;
        'wX' - warp mesh by vector field X, scale by factor F;
        'fF[%%]' - set scale factor F for warp/glyphs, see --factor option;
        'cC' - select Cth field component;
        'r' - recalculate cell data to point data;
        'sT' - plot data in step T;
        'mM' - plot cells with mat_id=M;
        'vP' - set plotting style to P: s=surface, w=wireframe, p=points;
        'e' - print edges;
        'oV' - set opacity to V;
        'pI' - plot in slot I""",
    'fields_map':
        'map fields and cell groups, e.g. 1:u1,p1 2:u2,p2',
    'outline':
        'plot mesh outline',
    'warp':
        'warp mesh by vector field',
    'factor':
        'scaling factor F for mesh warp and glyphs.'
        ' Append "%%" to scale relatively to the minimum bounding box size.',
    'edges':
        'plot cell edges',
    'isosurfaces':
        'plot isosurfaces [default: %(default)s]',
    'opacity':
        'set opacity [default: %(default)s]',
    'color_map':
        'set color_map, e.g. hot, cool, bone, etc. [default: %(default)s]',
    'color_limits':
        'set color bar range (min, max) for scalars'
        ' [default: given by data limits]',
    'axes_options':
        'options for directional axes, e.g. xlabel="z1" ylabel="z2",'
        ' zlabel="z3"',
    'no_axes':
        'hide orientation axes',
    'no_scalar_bars':
        'hide scalar bars',
    'grid_vector1':
        'define positions of plots along grid axis 1 [default: "0, 0, 1.6"]',
    'grid_vector2':
        'define positions of plots along grid axis 2 [default: "0, 1.6, 0"]',
    'max_plots':
        'maximum number of plots along grid axis 1'
        ' [default: 4]',
    'view':
        'camera azimuth, elevation angles, and optionally zoom factor'
        ' [default: "225,75,0.9"]',
    'camera_position':
        'define camera position',
    'window_size':
        'define size of plotting window',
    'animation':
        'create animation, mp4 file type supported',
    'framerate':
        'set framerate for animation',
    'screenshot':
        'save screenshot to file',
    'off_screen':
        'off screen plots, e.g. when screenshotting',
    'no_labels':
        'hide plot labels',
    'label_position':
        'define position of plot labels [default: "-1, -1, 0, 0.2"]',
    'scalar_bar_size':
        'define size of scalar bars [default: "0.15, 0.05"]',
    'scalar_bar_position':
        'define position of scalar bars [default: "0.8, 0.02, 0, 1.5"]',
    'step':
        'select data in a given time step',
    '2d_view':
        '2d view of XY plane',
    'probes':
        'visualize probes in the given files',
    'no_probe_labels':
        'hide probe labels',
    'no_step_time':
        'hide current step and time',
}


def main():
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('-f', '--fields', metavar='field_spec',
                        action=FieldOptsToListAction, nargs="+", dest='fields',
                        default=[], help=helps['fields'])
    parser.add_argument('--fields-map', metavar='map',
                        action=FieldOptsToListAction, nargs="+",
                        dest='fields_map',
                        default=[], help=helps['fields_map'])
    parser.add_argument('-s', '--step', metavar='step',
                        action=StoreNumberAction, dest='step',
                        default=0, help=helps['step'])
    parser.add_argument('-l', '--outline',
                        action='store_true', dest='outline',
                        default=False, help=helps['outline'])
    parser.add_argument('-i', '--isosurfaces',
                        action='store', dest='isosurfaces',
                        default=0, help=helps['isosurfaces'])
    parser.add_argument('-e', '--edges',
                        action='store_true', dest='show_edges',
                        default=False, help=helps['edges'])
    parser.add_argument('-w', '--warp', metavar='field',
                        action='store', dest='warp',
                        default=None, help=helps['warp'])
    parser.add_argument('--factor', metavar='factor',
                        action=StoreNumberAction, dest='factor',
                        default=1., help=helps['factor'])
    parser.add_argument('--opacity', metavar='opacity',
                        action=StoreNumberAction, dest='opacity',
                        default=1., help=helps['opacity'])
    parser.add_argument('--color-map', metavar='cmap',
                        action='store', dest='color_map',
                        default='viridis', help=helps['color_map'])
    parser.add_argument('--color-limits', metavar='min,max',
                        action=StoreNumberAction, dest='color_limits',
                        default=None, help=helps['color_limits'])
    parser.add_argument('--axes-options', metavar='options',
                        action=OptsToListAction, nargs="+",
                        dest='axes_options',
                        default=[], help=helps['axes_options'])
    parser.add_argument('--no-axes',
                        action='store_false', dest='axes_visibility',
                        default=True, help=helps['no_axes'])
    parser.add_argument('--grid-vector1', metavar='grid_vector1',
                        action=StoreNumberAction, dest='grid_vector1',
                        default=None, help=helps['grid_vector1'])
    parser.add_argument('--grid-vector2', metavar='grid_vector2',
                        action=StoreNumberAction, dest='grid_vector2',
                        default=None, help=helps['grid_vector2'])
    parser.add_argument('--max-plots',
                        action=StoreNumberAction, dest='max_plots',
                        default=4, help=helps['max_plots'])
    parser.add_argument('--no-labels',
                        action='store_false', dest='show_labels',
                        default=True, help=helps['no_labels'])
    parser.add_argument('--label-position', metavar='position',
                        action=StoreNumberAction, dest='label_position',
                        default=[-1, -1, 0, 0.2], help=helps['label_position'])
    parser.add_argument('--no-scalar-bars',
                        action='store_false', dest='show_scalar_bars',
                        default=True, help=helps['no_scalar_bars'])
    parser.add_argument('--scalar-bar-size', metavar='size',
                        action=StoreNumberAction, dest='scalar_bar_size',
                        default=[0.15, 0.05],
                        help=helps['scalar_bar_size'])
    parser.add_argument('--scalar-bar-position', metavar='position',
                        action=StoreNumberAction, dest='scalar_bar_position',
                        default=[0.8, 0.02, 0, 1.5],
                        help=helps['scalar_bar_position'])
    parser.add_argument('-v', '--view', metavar='position',
                        action=StoreNumberAction, dest='camera',
                        default=None, help=helps['view'])
    parser.add_argument('--camera-position', metavar='camera_position',
                        action=StoreNumberAction, dest='camera_position',
                        default=None, help=helps['camera_position'])
    parser.add_argument('--window-size', metavar='window_size',
                        action=StoreNumberAction, dest='window_size',
                        default=pv.global_theme.window_size,
                        help=helps['window_size'])
    parser.add_argument('-a', '--animation', metavar='output_file',
                        action='store', dest='anim_output_file',
                        default=None, help=helps['animation'])
    parser.add_argument('-r', '--frame-rate', metavar='rate',
                        action=StoreNumberAction, dest='framerate',
                        default=2.5, help=helps['framerate'])
    parser.add_argument('-o', '--screenshot', metavar='output_file',
                        action='store', dest='screenshot',
                        default=None, help=helps['screenshot'])
    parser.add_argument('--off-screen',
                        action='store_true', dest='off_screen',
                        default=False, help=helps['off_screen'])
    parser.add_argument('-2', '--2d-view',
                        action='store_true', dest='view_2d',
                        default=False, help=helps['2d_view'])
    parser.add_argument('--probes', metavar='filename', nargs='+',
                        dest='probes', default=None, help=helps['probes'])
    parser.add_argument('--no-probe-labels',
                        action='store_false', dest='show_probe_labels',
                        default=True, help=helps['no_probe_labels'])
    parser.add_argument('--no-step-time',
                        action='store_false', dest='show_step_time',
                        default=True, help=helps['no_step_time'])

    parser.add_argument('filenames', nargs='+')
    options = parser.parse_args()

    if options.probes is not None:
        annotations = read_probes_as_annotations(options.probes,
                                                 options.show_probe_labels)

    else:
        annotations = None

    anim_filename = options.anim_output_file
    if anim_filename:
        if anim_filename.endswith('.png') or anim_filename.endswith('.mp4'):
            options.off_screen = True

    pv.set_plot_theme("document")
    plotter = pv.Plotter(off_screen=options.off_screen,
                         title=make_title(options.filenames))

    if options.anim_output_file:
        _, _, n_steps = read_mesh(options.filenames, ret_n_steps=True)
        # dry run
        scalar_bar_limits = None
        if options.axes_visibility:
            plotter.add_axes(**dict(options.axes_options))
        for step in range(n_steps):
            plotter.clear()
            plotter, sb_limits = pv_plot(options.filenames, options,
                                         plotter=plotter, step=step,
                                         annotations=annotations,
                                         ret_scalar_bar_limits=True)
            if scalar_bar_limits is None:
                scalar_bar_limits = {k: [] for k in sb_limits.keys()}

            for k, v in sb_limits.items():
                scalar_bar_limits[k].append(v)

        cpos = _get_cpos(plotter, options)

        if cpos is not None:
            plotter.camera_position = cpos

        elif options.view_2d:
            plotter.view_xy()

        if anim_filename.endswith('.png'):
            from sfepy.base.ioutils import edit_filename
            fig_name = edit_filename(anim_filename, suffix='{step:05d}')
            plotter.show(auto_close=False)

        elif anim_filename.endswith('.gif'):
            from sfepy.base.ioutils import edit_filename
            fig_name = None
            plotter.open_gif(anim_filename)

        else:
            fig_name = None
            plotter.open_movie(anim_filename, options.framerate)
            plotter.show(auto_close=False)

        for k in scalar_bar_limits.keys():
            lims = scalar_bar_limits[k]
            clim = (nm.min([v[0] for v in lims]),
                    nm.max([v[1] for v in lims]))
            scalar_bar_limits[k] = clim

        # plot frames
        for step in range(n_steps):
            plotter.clear()
            plotter = pv_plot(options.filenames, options, plotter=plotter,
                              step=step, annotations=annotations,
                              scalar_bar_limits=scalar_bar_limits)
            if options.axes_visibility:
                plotter.add_axes(**dict(options.axes_options))

            if fig_name is None:
                plotter.write_frame()

            else:
                plotter.screenshot(fig_name.format(step=step), return_img=False)

        plotter.close()
    else:
        plotter = pv_plot(options.filenames, options, plotter=plotter,
                          annotations=annotations)
        if options.axes_visibility:
            plotter.add_axes(**dict(options.axes_options))

        plotter.add_key_event(
            'Prior', lambda: pv_plot(options.filenames,
                                     options,
                                     step=plotter.resview_step,
                                     step_inc=-1,
                                     annotations=annotations,
                                     plotter=plotter)
        )
        plotter.add_key_event(
            'Next', lambda: pv_plot(options.filenames,
                                    options,
                                    step=plotter.resview_step,
                                    step_inc=1,
                                    annotations=annotations,
                                    plotter=plotter)
        )

        # Does not work for meshes with no z component.
        plotter.add_key_event(
            'c', lambda: print_camera_position(plotter)
        )

        cpos = _get_cpos(plotter, options)
        if (cpos is None) and options.view_2d:
            plotter.view_xy()

        plotter.show(cpos=cpos, screenshot=options.screenshot,
                     window_size=options.window_size)

        if options.screenshot is not None and osp.exists(options.screenshot):
            print(f'saved: {options.screenshot}')

if __name__ == '__main__':
    main()

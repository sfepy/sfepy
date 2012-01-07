import os
import shutil, tempfile

import numpy as nm

try:
    from enthought.traits.api \
         import HasTraits, Instance, Button, Int, Bool, on_trait_change
    from enthought.traits.ui.api \
         import View, Item, Group, HGroup, spring
    from  enthought.traits.ui.editors.range_editor import RangeEditor
    from enthought.tvtk.pyface.scene_editor import SceneEditor
    from enthought.mayavi.tools.mlab_scene_model import MlabSceneModel
    from enthought.mayavi.core.ui.mayavi_scene import MayaviScene

except ImportError:
    from traits.api \
         import HasTraits, Instance, Button, Int, Bool, on_trait_change
    from traitsui.api \
         import View, Item, Group, HGroup, spring
    from traitsui.editors.range_editor import RangeEditor
    from tvtk.pyface.scene_editor import SceneEditor
    from mayavi.tools.mlab_scene_model import MlabSceneModel
    from mayavi.core.ui.mayavi_scene import MayaviScene


from dataset_manager import DatasetManager

from sfepy.base.base import (insert_as_static_method, output,
                             get_arguments, get_default, Struct, basestr)
from sfepy.linalg import cycle
from sfepy.solvers.ts import get_print_info
from sfepy.postprocess.utils import mlab
from sfepy.postprocess.sources import create_file_source, FileSource

def get_glyphs_scale_factor(rng, rel_scaling, bbox):
    delta = rng[1] - rng[0]
    dx = nm.max((bbox[1,:] - bbox[0,:]))
    if rel_scaling is None:
        rel_scaling = 0.02 # -> delta fits 50x into dx.
    return rel_scaling * dx / delta

def add_surf(obj, position, opacity=1.0):
    surf = mlab.pipeline.surface(obj, opacity=opacity)
    surf.actor.actor.position = position
    return surf

def add_scalar_cut_plane(obj, position, normal, opacity=1.0):
    scp = mlab.pipeline.scalar_cut_plane(obj, opacity=opacity)
    scp.actor.actor.position = position
    scp.implicit_plane.visible = False
    scp.implicit_plane.normal = normal

    return scp

def add_vector_cut_plane(obj, position, normal, bbox, rel_scaling=None,
                         scale_factor='auto', clamping=False,
                         opacity=1.0):
    vcp = mlab.pipeline.vector_cut_plane(obj, opacity=opacity)
    if scale_factor == 'auto':
        rng = vcp.glyph.glyph.range
        scale_factor = get_glyphs_scale_factor(rng, rel_scaling, bbox)

    vcp.glyph.color_mode = 'color_by_vector'
    vcp.glyph.scale_mode = 'scale_by_vector'
    vcp.glyph.glyph.clamping = clamping
    vcp.glyph.glyph.scale_factor = scale_factor
    vcp.glyph.glyph_source.glyph_position = 'tail'

    vcp.actor.actor.position = position
    vcp.implicit_plane.normal = normal
    vcp.implicit_plane.widget.enabled = False

    return vcp

def add_iso_surface(obj, position, contours=10, opacity=1.0):
    obj = mlab.pipeline.iso_surface(obj, contours=contours, opacity=opacity)
    obj.actor.actor.position = position
    return obj

def add_subdomains_surface(obj, position, mat_id_name='mat_id',
                           threshold_limits=(None, None), **kwargs):
    dm = DatasetManager(dataset=obj.outputs[0])
    mat_id = dm.cell_scalars[mat_id_name]

    rm = mat_id.min(), mat_id.max()

    active = mlab.pipeline.set_active_attribute(obj)
    active.cell_scalars_name = mat_id_name

    aa = mlab.pipeline.set_active_attribute(obj)
    aa.cell_scalars_name = mat_id_name

    threshold = mlab.pipeline.threshold(aa)
    threshold.threshold_filter.progress = 1.0
    if threshold_limits[0] is not None:
        threshold.lower_threshold = threshold_limits[0] + 0.1
    if threshold_limits[1] is not None:
        threshold.upper_threshold = threshold_limits[1] - 0.1

    threshold.auto_reset_lower = False
    threshold.auto_reset_upper = False

    surface = mlab.pipeline.surface(threshold, opacity=0.3)
    surface.actor.actor.position = position

    module_manager = surface.parent
    lm = module_manager.scalar_lut_manager
    lm.lut_mode = 'Blues'
    if (rm[1] - rm[0]) == 1:
        lm.reverse_lut = True

    surface2 = mlab.pipeline.surface(active, opacity=0.2)
    surface2.actor.actor.position = position

    module_manager = surface2.parent
    module_manager.scalar_lut_manager.lut_mode = 'Blues'

    return surface, surface2

def add_glyphs(obj, position, bbox, rel_scaling=None,
               scale_factor='auto', clamping=False, color=None):

    glyphs = mlab.pipeline.glyph(obj, mode='2darrow', scale_mode='vector',
                                 color=color, opacity=1.0) 
    if scale_factor == 'auto':
        rng = glyphs.glyph.glyph.range
        scale_factor = get_glyphs_scale_factor(rng, rel_scaling, bbox)

    glyphs.glyph.color_mode = 'color_by_vector'
    glyphs.glyph.scale_mode = 'scale_by_vector'
    glyphs.glyph.glyph.clamping = clamping
    glyphs.glyph.glyph.scale_factor = scale_factor
    glyphs.glyph.glyph_source.glyph_position = 'tail'
    glyphs.actor.actor.position = position
    return glyphs

def add_text(obj, position, text, width=None, color=(0, 0, 0)):
    if width is None:
        width = 0.02 * len(text)
    t = mlab.text(x=position[0], y=position[1], text=text,
                  z=position[2], color=color, width=width)
    return t

def get_position_counts(n_data, layout):
    n_col = max(1.0, min(5.0, nm.fix(nm.sqrt(n_data))))
    n_row = int(nm.ceil(n_data / n_col))
    n_col = int(n_col)
    if layout == 'rowcol':
        n_row, n_col = n_col, n_row
    elif layout == 'row':
        n_row, n_col = 1, n_data
    elif layout == 'col':
        n_row, n_col = n_data, 1
    else: # layout == 'rowcol':
        pass
    return n_row, n_col

class Viewer(Struct):
    """Class to automate visualization of various data using Mayavi. It can be
    used via postproc.py or isfepy the most easily.

    It can use any format that mlab.pipeline.open() handles, e.g. a VTK format.
    After opening a data file, all data (point, cell, scalars, vectors,
    tensors) are plotted in a grid layout.

    Parameters:

    watch : bool
        If True, watch the file for changes and update the mayavi
        pipeline automatically.
    animate : bool
        If True, save a view snaphost for each time step and exit.
    anim_format : str
        If set to a ffmpeg-supported format (e.g. mov, avi, mpg), ffmpeg is
        installed and results of multiple time steps are given, an animation is
        created in the same directory as the view images.
    ffmpeg_options : str
        The ffmpeg animation encoding options.
    output_dir : str
        The output directory, where view snapshots will be saved.

    Examples:

    >>> view = Viewer('file.vtk')
    >>> view() # view with default parameters
    >>> view(layout='col') # use column layout
    """
    def __init__(self, filename, watch=False,
                 animate=False, anim_format=None, ffmpeg_options=None,
                 output_dir='.', offscreen=False, auto_screenshot=True):
        Struct.__init__(self,
                        filename = filename,
                        watch = watch,
                        animate = animate,
                        anim_format = anim_format,
                        ffmpeg_options = ffmpeg_options,
                        output_dir = output_dir,
                        offscreen = offscreen,
                        auto_screenshot = auto_screenshot,
                        scene = None,
                        gui = None)
        self.options = get_arguments(omit = ['self'])

        if mlab is None:
            output('mlab cannot be imported, check your installation!')
            insert_as_static_method(self.__class__, '__call__', self.call_empty)
        else:
            insert_as_static_method(self.__class__, '__call__', self.call_mlab)

    def get_data_names(self, source=None, detailed=False):
        if source is None:
            mlab.options.offscreen = self.offscreen
            mlab.figure(fgcolor=self.fgcolor, bgcolor=self.bgcolor,
                        size=(1, 1))
            source = mlab.pipeline.open(self.filename)
        point_scalar_names = sorted( source._point_scalars_list[:-1] )
        point_vector_names = sorted( source._point_vectors_list[:-1] )
        point_tensor_names = sorted( source._point_tensors_list[:-1] )
        cell_scalar_names = sorted( source._cell_scalars_list[:-1] )
        cell_vector_names = sorted( source._cell_vectors_list[:-1] )
        cell_tensor_names = sorted( source._cell_tensors_list[:-1] )

        p_names = [['point', 'scalars', name] for name in point_scalar_names]
        p_names += [['point', 'vectors', name] for name in point_vector_names]
        p_names += [['point', 'tensors', name] for name in point_tensor_names]
        c_names = [['cell', 'scalars', name] for name in cell_scalar_names]
        c_names += [['cell', 'vectors', name] for name in cell_vector_names]
        c_names += [['cell', 'tensors', name] for name in cell_tensor_names]

        if detailed:
            return p_names, c_names
        else:
            return p_names + c_names

    def set_source_filename(self, filename):
        self.filename = filename
        try:
            self.file_source.set_filename(filename, self.scene.children[0])
        except (AttributeError, IndexError): # No sources yet.
            pass

    def save_image(self, filename):
        """Save a snapshot of the current scene."""
        if self.output_dir and not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        name = os.path.join(self.output_dir, filename)
        output('saving %s...' % name)
        self.scene.scene.save(name)
        output('...done')

    def get_animation_info(self, filename, add_output_dir=True, rng=None):
        if rng is None:
            rng = self.file_source.get_step_range()

        base, ext = os.path.splitext(filename)
        if add_output_dir:
            base = os.path.join(self.output_dir, base)

        n_digit, fmt, suffix = get_print_info(rng[1] - rng[0] + 1)
        return base, suffix, ext

    def save_animation(self, filename):
        """Animate the current scene view for all the time steps and save
        a snapshot of each step view."""
        rng = self.file_source.get_step_range()
        base, suffix, ext = self.get_animation_info(filename,
                                                    add_output_dir=False,
                                                    rng=rng)

        for step in xrange(rng[0], rng[1]+1):
            name = '.'.join((base, suffix % step, ext[1:]))
            output('%d: %s' % (step, name))

            self.set_step.step = step
            self.save_image(name)

    def encode_animation(self, filename, format, ffmpeg_options=None):
        if ffmpeg_options is None:
            ffmpeg_options = '-r 10 -sameq'

        base, suffix, ext = self.get_animation_info(filename)
        anim_name = '.'.join((base, format))
        cmd = 'ffmpeg %s -i %s %s' % (ffmpeg_options,
                                      '.'.join((base, suffix, ext[1:])),
                                      anim_name)
        output('creating animation "%s"...' % anim_name)
        try:
            os.system(cmd) 
        except:
            output('...warning: animation not created, is ffmpeg installed?')
        else:
            output('...done')

        return anim_name

    def get_size_hint(self, layout, resolution=None):
        if resolution is not None:
            size = resolution
        elif layout == 'rowcol':
            size = (800, 600)
        elif layout == 'row':
            size = (1000, 600)
        elif layout == 'col':
            size = (600, 1000)
        else:
            size = (600, 800)
        return size

    def build_mlab_pipeline(self, file_source=None, is_3d=False,
                            layout='rowcol',
                            scalar_mode='iso_surface',
                            vector_mode='arrows_norm',
                            rel_scaling=None, clamping=False,
                            ranges=None, is_scalar_bar=False,
                            is_wireframe=False, opacity=1.0,
                            subdomains_args=None,
                            rel_text_width=None,
                            filter_names=None, group_names=None,
                            only_names=None,
                            domain_specific=None, **kwargs):
        """Sets self.source, self.is_3d_data """
        file_source = get_default(file_source, self.file_source,
                                  'file_source not set!')
        filter_names = get_default(filter_names, [])
        domain_specific = get_default(domain_specific, {})

        if subdomains_args is not None:
            is_subdomains = True
            file_source.setup_mat_id(subdomains_args['mat_id_name'],
                                     subdomains_args['single_color'])

        else:
            is_subdomains = False

        self.source = source = self.file_source()

        # Hack to prevent mayavi switching to point scalar on source
        # change.
        if len(source._point_scalars_list):
            source.point_scalars_name = ''

        bbox = file_source.get_bounding_box()
        dx = 1.1 * (bbox[1,:] - bbox[0,:])

        float_eps = nm.finfo(nm.float64).eps
        self.is_3d_data = abs(dx[2]) > (10.0 * float_eps)

        p_names, c_names = self.get_data_names(source, detailed=True)
        names = p_names + c_names
        if only_names is None:
            names = [ii for ii in names if ii[2] not in filter_names]
        else:
            _names = [ii for ii in names if ii[2] in only_names]
            if len(_names) != len(only_names):
                output('warning: some names were not found!')
            if not len(_names):
                raise ValueError('no names were found! (%s not in %s)'
                                 % (only_names, [name[2] for name in names]))
            names = _names

        if group_names is not None:
            ndict = {}
            for name in names:
                ndict[name[2]] = name

            repeat = []
            _names = []
            aux = set(name[2] for name in names)
            for group in group_names:
                aux.difference_update(group)
                repeat.append(len(group))
                for name in group:
                    _names.append(ndict[name])

            repeat.extend([1] * len(aux))
            n_pos = len(repeat)

            names = _names
            n_data = len(names)

        else:
            n_pos = n_data = len(names)
            repeat = [1] * n_data

        def _make_iterator(repeat, n_row, n_col):
            ii = 0
            for ij, iric in enumerate(cycle((n_row, n_col))):
                ir, ic = iric
                if ij < len(repeat):
                    for ik in xrange(repeat[ij]):
                        yield ii, ir, ic
                        ii += 1

        n_row, n_col = get_position_counts(n_pos, layout)
        if layout[:3] == 'col':
            iterator = _make_iterator(repeat, n_col, n_row)
        else:
            iterator = _make_iterator(repeat, n_row, n_col)

        max_label_width = nm.max([len(ii[2]) for ii in names] + [5]) + 2

        if c_names:
            ctp = mlab.pipeline.cell_to_point_data(source)

        else:
            ctp = None

        self.scalar_bars = []

        for ii, ir, ic in iterator:
            if layout[:3] == 'col':
                ir, ic = ic, ir
            if ii == n_data: break
            family, kind, name = names[ii]

            is_magnitude = False
            position = nm.array([dx[0] * ic, dx[1] * (n_row - ir - 1), 0])
            output(family, kind, name, position)

            if name in domain_specific:
                ds = domain_specific[name]
                out = ds(source, ctp, bbox, position, family, kind, name)
                if len(out) == 4:
                    kind, name, active, bars = out
                    self.scalar_bars.extend(bars)

                else:
                    kind, name, active = out

            elif kind == 'scalars':
                if family == 'point':
                    active = mlab.pipeline.set_active_attribute(source)

                else:
                    if is_3d and ('iso_surface' in scalar_mode):
                        active = mlab.pipeline.set_active_attribute(ctp)

                    else:
                        active = mlab.pipeline.set_active_attribute(source)

                setattr(active, '%s_%s_name' % (family, kind), name)

                if is_3d:
                    if 'cut_plane' in scalar_mode:
                        scp = add_scalar_cut_plane(active,
                                                   position, [1, 0, 0],
                                                   opacity=0.5*opacity)
                        scp = add_scalar_cut_plane(active,
                                                   position, [0, 1, 0],
                                                   opacity=0.5*opacity)
                        scp = add_scalar_cut_plane(active,
                                                   position, [0, 0, 1],
                                                   opacity=0.5*opacity)
                    if 'iso_surface' in scalar_mode:
                        active.point_scalars_name = name
                        iso = add_iso_surface(active, position,
                                              opacity=0.3*opacity)
                else:
                    surf = add_surf(active, position, opacity=opacity)
                
            elif kind == 'vectors':
                if family == 'point':
                    active = mlab.pipeline.set_active_attribute(source)
                else:
                    active = mlab.pipeline.set_active_attribute(ctp)
                active.point_vectors_name = name

                if (ranges is not None) and (name in ranges):
                    sf = get_glyphs_scale_factor(ranges[name],
                                                 rel_scaling, bbox)
                else:
                    sf = None

                if 'arrows' in vector_mode:
                    glyphs = add_glyphs(active, position, bbox,
                                        rel_scaling=rel_scaling,
                                        clamping=clamping)
                    if sf is not None:
                        glyphs.glyph.glyph.scale_factor = sf
                        
                if 'warp' in vector_mode:
                    active = mlab.pipeline.warp_vector(active)
                    active.filter.scale_factor = rel_scaling

                if 'norm' in vector_mode:
                    active = mlab.pipeline.extract_vector_norm(active)
                    if 'arrows' in vector_mode:
                        a_opacity = 0.3 * opacity
                    else:
                        a_opacity = opacity
                    surf = add_surf(active, position, opacity=a_opacity)

                if is_3d:
                    if 'cut_plane' in vector_mode:
                        for normal in [[1, 0, 0], [0, 1, 0], [0, 0, 1]]:
                            vcp = add_vector_cut_plane(active,
                                                       position, normal, bbox,
                                                       rel_scaling=rel_scaling,
                                                       clamping=clamping,
                                                       opacity=0.5*opacity)
                            if sf is not None:
                                vcp.glyph.glyph.scale_factor = sf

            elif kind == 'tensors':
                if family == 'point':
                    active = mlab.pipeline.set_active_attribute(source)
                else:
                    active = mlab.pipeline.set_active_attribute(ctp)
                active.point_tensors_name = name

                active = mlab.pipeline.extract_tensor_components(active)
                is_magnitude = True
                if is_3d:
                    if 'cut_plane' in scalar_mode:
                        scp = add_scalar_cut_plane(active,
                                                   position, [1, 0, 0],
                                                   opacity=0.5)
                        scp = add_scalar_cut_plane(active,
                                                   position, [0, 1, 0],
                                                   opacity=0.5 )
                        scp = add_scalar_cut_plane(active,
                                                   position, [0, 0, 1],
                                                   opacity=0.5 )
                    if 'iso_surface' in scalar_mode:
                        iso = add_iso_surface(active, position,
                                              opacity=0.3*opacity)
                else:
                    surf = add_surf(active, position, opacity=opacity)

            else:
                raise ValueError('bad kind! (%s)' % kind)

            if (ranges is not None) and (name in ranges):
                mm = active.children[0]
                if (kind == 'scalars') or (kind == 'tensors'):
                    lm = mm.scalar_lut_manager

                else: # kind == 'vectors': 
                    lm = mm.vector_lut_manager
                    
                lm.use_default_range = False
                lm.data_range = ranges[name]

            if is_subdomains:
                add_subdomains_surface(source, position, **subdomains_args)

            if is_wireframe:
                surf = add_surf(source, position, opacity=opacity)
                surf.actor.property.representation = 'wireframe'
                surf.actor.mapper.scalar_visibility = False

            if is_scalar_bar:
                mm = active.children[0]
                if (kind == 'scalars') or (kind == 'tensors'):
                    lm = mm.scalar_lut_manager

                else: # kind == 'vectors': 
                    lm = mm.vector_lut_manager

                self.scalar_bars.append((family, name, lm))
                
            if rel_text_width > (10 * float_eps):
                position[2] = 0.5 * dx[2]
                if is_magnitude:
                    name = '|%s|' % name
                text = add_text(active, position, name,
                                float(rel_text_width * len(name))
                                / float(max_label_width), color=self.fgcolor)

        if not names:
            # No data, so just show the mesh.
            surf = add_surf(source, (0.0, 0.0, 0.0), opacity=opacity)
            surf.actor.property.color = (0.8, 0.8, 0.8)

            if is_subdomains:
                add_subdomains_surface(source, (0.0, 0.0, 0.0),
                                       **subdomains_args)

            if is_wireframe:
                surf = add_surf(source, (0.0, 0.0, 0.0), opacity=opacity)
                surf.actor.property.representation = 'wireframe'
                surf.actor.mapper.scalar_visibility = False

    def render_scene(self, scene, options):
        """
        Render the scene, preferably after it has been activated.
        """
        scene.scene.disable_render = True

        self.build_mlab_pipeline(**options)
        self.source.update() # Force source update to see e.g. streamlines.
        scene.scene.reset_zoom()

        view = options['view']
        if view is None:
            if options['is_3d'] or self.is_3d_data:
                self.view = (45, 45)
            else:
                self.view = (0, 0)
        else:
            self.view = view
        self.roll = options['roll']

        scene.scene.disable_render = False

        anti_aliasing = options['anti_aliasing']
        if anti_aliasing is not None:
            scene.scene.anti_aliasing_frames = anti_aliasing

    def show_scalar_bars(self, scalar_bars):
        for ii, (family, name, lm) in enumerate(scalar_bars):
            x0, y0 = 0.03, 1.0 - 0.01 - float(ii) / 15.0 - 0.07
            x1, y1 = 0.4, 0.07
            lm.scalar_bar_representation.position = [x0, y0]
            lm.scalar_bar_representation.position2 = [x1, y1]
            lm.number_of_labels = 5
            lm.scalar_bar.orientation = 'horizontal'
            lm.data_name = '%s: %s' % (family, name)
            lm.scalar_bar.title_text_property.font_size = 20
            lm.scalar_bar.label_text_property.font_size = 16
            lm.show_scalar_bar = True

    def reset_view(self):
        mlab.view(*self.view)
        mlab.roll(self.roll)
        self.scene.scene.camera.zoom(1.0)

    def __call__(self, *args, **kwargs):
        """
        This is either call_mlab() or call_empty().
        """
        pass
            
    def call_empty(self, *args, **kwargs):
        pass

    def call_mlab(self, scene=None, show=True, is_3d=False,
                  view=None, roll=None,
                  fgcolor=(1.0, 1.0, 1.0), bgcolor=(0.0, 0.0, 0.0),
                  layout='rowcol', scalar_mode='iso_surface',
                  vector_mode='arrows_norm', rel_scaling=None, clamping=False,
                  ranges=None, is_scalar_bar=False, is_wireframe=False,
                  opacity=1.0, subdomains_args=None, rel_text_width=None,
                  fig_filename='view.png', resolution = None,
                  filter_names=None, only_names=None, group_names=None, step=0,
                  anti_aliasing=None, domain_specific=None):
        """
        By default, all data (point, cell, scalars, vectors, tensors)
        are plotted in a grid layout, except data named 'node_groups',
        'mat_id' which are usually not interesting.

        Parameters
        ----------
        show : bool
            Call mlab.show().
        is_3d : bool
            If True, use scalar cut planes instead of surface for certain
            datasets. Also sets 3D view mode.
        view : tuple
            Azimuth, elevation angles, distance and focal point as in
            `mlab.view()`.
        roll : float
            Roll angle tuple as in mlab.roll().
        fgcolor : tuple of floats (R, G, B)
            The foreground color, that is the color of all text
            annotation labels (axes, orientation axes, scalar bar
            labels).
        bgcolor : tuple of floats (R, G, B)
            The background color.
        layout : str
            Grid layout for placing the datasets. Possible values are:
            'row', 'col', 'rowcol', 'colrow'.
        scalar_mode : str
             Mode for plotting scalars and tensor magnitudes, one of
             'cut_plane', 'iso_surface', 'both'. 
        vector_mode : str
             Mode for plotting vectors, one of 'arrows', 'norm', 'arrows_norm',
             'warp_norm'.
        rel_scaling : float
            Relative scaling of glyphs for vector datasets.
        clamping : bool
            Clamping for vector datasets.
        ranges : dict
            List of data ranges in the form {name : (min, max), ...}.
        is_scalar_bar : bool
            If True, show a scalar bar for each data.
        is_wireframe : bool
            If True, show a wireframe of mesh surface bar for each data.
        opacity : float
            Global surface and wireframe opacity setting in [0.0, 1.0],
        subdomains_args : tuple
            Tuple of (mat_id_name, threshold_limits, single_color), see
            :func:`add_subdomains_surface`, or None.
        rel_text_width : float
            Relative text width.
        fig_filename : str
            File name for saving the resulting scene figure, if
            self.auto_screenshot is True.
        resolution : tuple
            Scene and figure resolution. If None, it is set
            automatically according to the layout.
        filter_names : list of strings
            Omit the listed datasets. If None, it is initialized to
            ['node_groups', 'mat_id']. Pass [] if you need no filtering.
        only_names : list of strings
            Draw only the listed datasets. If None, it is initialized all names
            besides those in filter_names.
        group_names : list of tuples
            List of data names in the form [(name1, ..., nameN), (...)]. Plots
            of data named in each group are superimposed. Repetitions of names
            are possible.
        step : int
            The time step to display.
        anti_aliasing : int
            Value of anti-aliasing.
        domain_specific : dict
            Domain-specific drawing functions and configurations.
        """
        self.fgcolor = fgcolor
        self.bgcolor = bgcolor

        if filter_names is None:
            filter_names = ['node_groups', 'mat_id']

        if rel_text_width is None:
            rel_text_width = 0.02

        if isinstance(scalar_mode, basestr):
            if scalar_mode == 'both':
                scalar_mode = ('cut_plane', 'iso_surface')
            elif scalar_mode in ('cut_plane', 'iso_surface'):
                scalar_mode = (scalar_mode,)
            else:
                raise ValueError('bad value of scalar_mode parameter! (%s)'
                                 % scalar_mode)
        else:
            for sm in scalar_mode:
                if not sm in ('cut_plane', 'iso_surface'):
                    raise ValueError('bad value of scalar_mode parameter! (%s)'
                                     % sm)

        if isinstance(vector_mode, basestr):
            if vector_mode == 'arrows_norm':
                vector_mode = ('arrows', 'norm')
            elif vector_mode == 'warp_norm':
                vector_mode = ('warp', 'norm')
            elif vector_mode in ('arrows', 'norm'):
                vector_mode = (vector_mode,)
            elif vector_mode == 'cut_plane':
                if is_3d:
                    vector_mode = ('cut_plane',)
                else:
                    vector_mode = ('arrows',)
            else:
                raise ValueError('bad value of vector_mode parameter! (%s)'
                                 % vector_mode)
        else:
            for vm in vector_mode:
                if not vm in ('arrows', 'norm', 'warp'):
                    raise ValueError('bad value of vector_mode parameter! (%s)'
                                     % vm)

        mlab.options.offscreen = self.offscreen

        self.size_hint = self.get_size_hint(layout, resolution=resolution)

        is_new_scene = False

        if scene is not None:
            if scene is not self.scene:
                is_new_scene = True
                self.scene = scene
            gui = None

        else:
            if (self.scene is not None) and (not self.scene.running):
                self.scene = None

            if self.scene is None:
                if self.offscreen:
                    gui = None
                    scene = mlab.figure(fgcolor=fgcolor, bgcolor=bgcolor,
                                        size=self.size_hint)

                else:
                    gui = ViewerGUI(viewer=self,
                                    fgcolor=fgcolor, bgcolor=bgcolor)
                    scene = gui.scene.mayavi_scene

                if scene is not self.scene:
                    is_new_scene = True
                    self.scene = scene

            else:
                gui = self.gui
                scene = self.scene

        self.engine = mlab.get_engine()
        self.engine.current_scene = self.scene

        self.gui = gui

        self.file_source = create_file_source(self.filename, watch=self.watch,
                                              offscreen=self.offscreen)
        has_several_steps = (nm.diff(self.file_source.get_step_range()) > 0)[0]

        if gui is not None:
            gui.has_several_steps = has_several_steps

        if has_several_steps:
            self.set_step = set_step = SetStep()
            set_step._viewer = self
            set_step._source = self.file_source
            set_step.step = step
            self.file_source.setup_notification(set_step, 'file_changed')

            if gui is not None:
                gui.set_step = set_step

        self.options.update(get_arguments(omit = ['self', 'file_source']))

        if gui is None:
            self.render_scene(scene, self.options)
            self.reset_view()
            if is_scalar_bar:
                self.show_scalar_bars(self.scalar_bars)

            if self.animate:
                self.save_animation(fig_filename)

            else:
                self.save_image(fig_filename)

        else:
            traits_view = View(
                Item('scene', editor=SceneEditor(scene_class=MayaviScene), 
                     show_label=False,
                     width=self.size_hint[0], height=self.size_hint[1],
                     style='custom',
                ),
                Group(Item('set_step', defined_when='set_step is not None',
                           show_label=False, style='custom'),
                ),
                HGroup(spring,
                       Item('button_view', show_label=False),
                       Item('button_make_snapshots', show_label=False,
                            enabled_when='has_several_steps == True'),
                       Item('button_make_animation', show_label=False,
                            enabled_when='has_several_steps == True'),),
                resizable=True,
                buttons=['OK'],
            )

            if is_new_scene:
                if show:
                    gui.configure_traits(view=traits_view)

                else:
                    gui.edit_traits(view=traits_view)

                    if self.auto_screenshot:
                        self.save_image(fig_filename)

        return gui

class SetStep(HasTraits):

    _viewer = Instance(Viewer)
    _source = Instance(FileSource)
    _step_editor = RangeEditor(low_name='step_low',
                               high_name='step_high',
                               label_width=28,
                               auto_set=True,
                               mode='slider')
    step = None
    step_low = Int
    step_high = Int
    file_changed = Bool(False)

    traits_view = View(
        Item('step', defined_when='step is not None',
             editor=_step_editor),
    )

    def __source_changed(self, old, new):
        rng = self._source.get_step_range()
        self.add_trait('step', Int(0))
        self.step_low, self.step_high = [int(ii) for ii in rng]

    def _step_changed(self, old, new):
        self._source.set_step(self.step)
        self._viewer.set_source_filename(self._source.filename)

    def _file_changed_changed(self, old, new):
        if new == True:
            rng = self._source.get_step_range()
            self.step_low, self.step_high = [int(ii) for ii in rng]

        self.file_changed = False

def make_animation(filename, view, roll, anim_format, options,
                   reuse_viewer=None):
    output_dir = tempfile.mkdtemp()

    viewer = Viewer(filename, watch=options.watch,
                    animate=True,
                    output_dir=output_dir,
                    offscreen=True)

    if reuse_viewer is None:
        viewer(show=False, is_3d=options.is_3d, view=view,
               roll=roll, layout=options.layout,
               scalar_mode=options.scalar_mode,
               vector_mode=options.vector_mode,
               rel_scaling=options.rel_scaling,
               clamping=options.clamping, ranges=options.ranges,
               is_scalar_bar=options.is_scalar_bar,
               is_wireframe=options.is_wireframe,
               opacity=options.opacity,
               subdomains_args=options.subdomains_args,
               rel_text_width=options.rel_text_width,
               fig_filename=options.fig_filename, resolution=options.resolution,
               filter_names=options.filter_names, only_names=options.only_names,
               group_names=options.group_names,
               anti_aliasing=options.anti_aliasing,
               domain_specific=options.domain_specific)

    else:
        viewer.file_source = reuse_viewer.file_source
        viewer.scene = reuse_viewer.scene
        viewer.set_step = reuse_viewer.set_step
        viewer.save_animation(options.fig_filename)

    op = os.path
    if anim_format != 'png':
        anim_name = viewer.encode_animation(options.fig_filename, anim_format,
                                            options.ffmpeg_options)

        shutil.move(anim_name, op.join(options.output_dir,
                                       op.split(anim_name)[1]))
        shutil.rmtree(output_dir)

    else:
        shutil.move(output_dir, op.join(options.output_dir, 'snapshots'))

    if reuse_viewer is None:
        mlab.close(viewer.scene)

class ViewerGUI(HasTraits):

    scene = Instance(MlabSceneModel, ())

    has_several_steps = Bool(False)
    viewer = Instance(Viewer)
    set_step = Instance(SetStep)
    button_view = Button('print view')
    button_make_animation = Button('make animation')
    button_make_snapshots = Button('make snapshots')

##     anim_process = Instance(Process)

    @on_trait_change('scene.activated')
    def _post_init(self, name, old, new):
        viewer = self.viewer
        viewer.render_scene(self.scene, viewer.options)
        viewer.reset_view()
        viewer.show_scalar_bars(viewer.scalar_bars)

    def __init__(self, fgcolor=(0.0, 0.0, 0.0), bgcolor=(1.0, 1.0, 1.0),
                 **traits):
        HasTraits.__init__(self, **traits)
        scene = self.scene.scene

        scene.foreground = fgcolor
        scene.background = bgcolor

    def _button_view_fired(self):
        self.scene.camera.print_traits()
        print 'view:', mlab.view()
        print 'roll:', mlab.roll()

    def _button_make_animation_fired(self):
        view = mlab.view()
        roll = mlab.roll()

##         if self.anim_process and self.anim_process.is_alive():
##             output('terminating previous animation process...')
##             self.anim_process.terminate()
##             output('...done')

##         output('starting animation process...')
##         self.anim_process = Process(target=make_animation,
##                                     args=(self.viewer.filename,
##                                           view[:2],
##                                           roll,
##                                           Struct(**self.viewer.options)))
##         self.anim_process.daemon = True
##         self.anim_process.start()
##         output('...done')

        make_animation(self.viewer.filename,
                       view,
                       roll,
                       'avi',
                       Struct(**self.viewer.options),
                       self.viewer)

    def _button_make_snapshots_fired(self):
        view = mlab.view()
        roll = mlab.roll()

        make_animation(self.viewer.filename,
                       view,
                       roll,
                       'png',
                       Struct(**self.viewer.options),
                       self.viewer)

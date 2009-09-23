from sfepy.base.base import *
from sfepy.base.la import cycle
from sfepy.postprocess.utils import mlab
from sfepy.postprocess.sources import create_file_source, FileSource

from enthought.traits.api \
     import HasTraits, Instance, Range, on_trait_change
from enthought.traits.ui.api \
     import View, Item, Group
from  enthought.traits.ui.editors.range_editor import RangeEditor
from enthought.tvtk.pyface.scene_editor import SceneEditor
from enthought.mayavi.tools.mlab_scene_model import MlabSceneModel
from enthought.mayavi.core.ui.mayavi_scene import MayaviScene

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

def add_iso_surface(obj, position, contours=10, opacity=1.0):
    obj = mlab.pipeline.iso_surface(obj, contours=contours, opacity=opacity)
    obj.actor.actor.position = position
    return obj
    
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
    n_col = min(5.0, nm.fix(nm.sqrt(n_data)))
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

    Examples
    --------
    >>> view = Viewer('file.vtk')
    >>> view() # view with default parameters
    >>> view(layout='col') # use column layout
    """
    def __init__(self, filename, output_dir='.', offscreen=False,
                 auto_screenshot=True):
        Struct.__init__(self,
                        filename = filename,
                        output_dir = output_dir,
                        offscreen = offscreen,
                        auto_screenshot = auto_screenshot,
                        mlab = mlab)

        if mlab is None:
            output('mlab cannot be imported, check your installation!')
            insert_as_static_method(self.__class__, '__call__', self.call_empty)
        else:
            insert_as_static_method(self.__class__, '__call__', self.call_mlab)
            
    def get_data_names(self, source=None, detailed=False):
        if source is None:
            mlab.options.offscreen = self.offscreen
            scene = mlab.figure(bgcolor=(1,1,1), fgcolor=(0, 0, 0), size=(1, 1))
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
        except IndexError: # No sources yet.
            pass

    def save_image(self, filename):
        name = os.path.join(self.output_dir, filename)
        output('saving %s...' % name)
        self.scene.scene.save(name)
        output('...done')

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

    def build_mlab_pipeline(self, file_source=None, is_3d=False, layout='rowcol',
                            scalar_mode='iso_surface',
                            rel_scaling=None, clamping=False,
                            ranges=None, is_scalar_bar=False,
                            rel_text_width=None,
                            filter_names=None, only_names=None):
        """Sets self.source, self.is_3d_data """
        file_source = get_default(file_source, self.file_source,
                                  'file_source not set!')

        if filter_names is None:
            filter_names = []

        self.source = source = self.file_source()

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
        n_data = len(names)
        n_row, n_col = get_position_counts(n_data, layout)

        max_label_width = nm.max([len(ii[2]) for ii in names] + [5]) + 2

        if c_names:
            ctp = mlab.pipeline.cell_to_point_data(source)

        if layout[:3] == 'col':
            iterator = enumerate(cycle((n_col, n_row)))
        else:
            iterator = enumerate(cycle((n_row, n_col)))

        self.scalar_bars = []

        for ii, (ir, ic) in iterator:
            if layout[:3] == 'col':
                ir, ic = ic, ir
            if ii == n_data: break
            family, kind, name = names[ii]

            is_magnitude = False
            position = nm.array([dx[0] * ic, dx[1] * (n_row - ir - 1), 0])
            
            output(family, kind, name, position)
            if kind == 'scalars':
                active = mlab.pipeline.set_active_attribute(source)
#                active.point_scalars_name = name
                setattr(active, '%s_%s_name' % (family, kind), name)

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
                        iso = add_iso_surface(active, position, opacity=0.3)
                else:
                    surf = add_surf(active, position)
                
            elif kind == 'vectors':
                if family == 'point':
                    active = mlab.pipeline.set_active_attribute(source)
                else:
                    active = mlab.pipeline.set_active_attribute(ctp)
                active.point_vectors_name = name

                glyphs = add_glyphs(active, position, bbox,
                                    rel_scaling=rel_scaling, clamping=clamping)

                if (ranges is not None) and (name in ranges):
                    sf = get_glyphs_scale_factor(ranges[name], rel_scaling, bbox)
                    glyphs.glyph.glyph.scale_factor = sf

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
                        iso = add_iso_surface(active, position, opacity=0.3)
                else:
                    surf = add_surf(active, position)

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
                                / float(max_label_width))

        if not names:
            # No data, so just show the mesh.
            surf = add_surf(source, (0.0, 0.0, 0.0))
            surf.actor.property.color = (0.8, 0.8, 0.8)

            surf = add_surf(source, (0.0, 0.0, 0.0))
            surf.actor.property.representation = 'wireframe'

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
        
    def __call__(self, *args, **kwargs):
        """
        This is either call_mlab() or call_empty().
        """
        pass
            
    def call_empty(self, *args, **kwargs):
        pass
    
    def call_mlab(self, scene=None, show=True, is_3d=False, view=None, roll=None,
                  layout='rowcol', scalar_mode='iso_surface',
                  rel_scaling=None, clamping=False,
                  ranges=None, is_scalar_bar=False, rel_text_width=None,
                  fig_filename='view.png', resolution = None,
                  filter_names=None, only_names=None, anti_aliasing=None):
        """By default, all data (point, cell, scalars, vectors, tensors) are
        plotted in a grid layout, except data named 'node_groups', 'mat_id' which
        are usually not interesting.

        Parameters
        ----------
        show : bool
            Call mlab.show().
        is_3d : bool
            If True, use scalar cut planes instead of surface for certain
            datasets. Also sets 3D view mode.
        view : tuple
            Azimuth, elevation angles as in mlab.view().
        roll : float
            Roll angle tuple as in mlab.roll().
        layout : str
            Grid layout for placing the datasets. Possible values are:
            'row', 'col', 'rowcol', 'colrow'.
        scalar_mode : str
             Mode for plotting scalars and tensor magnitudes, one of
             'cut_plane', 'iso_surface', 'both'. 
        rel_scaling : float
            Relative scaling of glyphs for vector datasets.
        clamping : bool
            Clamping for vector datasets.
        ranges : dict
            List of data ranges in the form {name : (min, max), ...}.
        is_scalar_bar : bool
            If True, show a scalar bar for each data.
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
        anti_aliasing : int
            Value of anti-aliasing.
        """
        if filter_names is None:
            filter_names = ['node_groups', 'mat_id']

        if rel_text_width is None:
            rel_text_width = 0.02

        if scalar_mode == 'both':
            scalar_mode = ('cut_plane', 'iso_surface')
        elif scalar_mode in ('cut_plane', 'iso_surface'):
            scalar_mode = (scalar_mode,)
        else:
            raise ValueError('bad value of scalar_mode parameter! (%s)'
                             % scalar_mode)

        mlab.options.offscreen = self.offscreen

        self.size_hint = self.get_size_hint(layout, resolution=resolution)

        if scene is None:
            if self.offscreen:
                self.scene = mlab.figure(bgcolor=(1,1,1), fgcolor=(0, 0, 0),
                                         size=self.size_hint)
                scene = self.scene
                gui = None

            else:
                gui = ViewerGUI(viewer=self)
                self.scene = scene = gui.scene.mlab.get_engine().current_scene

        else:
            gui = None
            self.scene = scene

        scene.scene.disable_render = True

        self.file_source = create_file_source(self.filename,
                                              offscreen=self.offscreen)
        if nm.diff(self.file_source.get_step_range()) > 0:
            set_step = SetStep()
            set_step._viewer = self
            set_step._source = self.file_source
            set_step.step = 0
            gui.set_step = set_step

        self.build_mlab_pipeline(is_3d=is_3d,
                                 layout=layout,
                                 scalar_mode=scalar_mode,
                                 rel_scaling=rel_scaling,
                                 clamping=clamping,
                                 ranges=ranges,
                                 is_scalar_bar=is_scalar_bar,
                                 rel_text_width=rel_text_width,
                                 filter_names=filter_names,
                                 only_names=only_names)
            
        scene.scene.reset_zoom()

        if view is None:
            if is_3d or self.is_3d_data:
                self.view = (45, 45)
            else:
                self.view = (0, 0)
        else:
            self.view = view
        self.roll = roll
      
        scene.scene.disable_render = False
        if anti_aliasing is not None:
            scene.scene.anti_aliasing_frames = anti_aliasing

## This works only after the scene is active!
##         if self.auto_screenshot:
##             self.save_image(fig_filename)

##         if is_scalar_bar:
##             self.show_scalar_bars(self.scalar_bars)

        if gui is None:
            mlab.view(*self.view)
            mlab.roll(self.roll)
            scene.scene.camera.zoom(1.0)
            if is_scalar_bar:
                self.show_scalar_bars(self.scalar_bars)
            self.save_image(fig_filename)

        else:
            traits_view = View(
                Item('scene', editor=SceneEditor(scene_class=MayaviScene), 
                     show_label=False,
                     width=self.size_hint[0], height=self.size_hint[1],
                     style='custom',
                ),
                Group(Item('set_step', defined_when='set_step is not None',
                           show_label=False, style='custom'),),
                resizable=True,
                buttons=['OK'],
            )

            if show:
                gui.configure_traits(view=traits_view)

            else:
                gui.edit_traits(view=traits_view)

        return gui

class SetStep(HasTraits):

    _viewer = Instance(Viewer)
    _source = Instance(FileSource)
    _step_editor = RangeEditor(low=0, high=0, auto_set=True, mode='slider')
    step = None

    traits_view = View(
        Item('step', defined_when='step is not None',
             editor=_step_editor),
    )

    def __source_changed(self, old, new):
        rng = self._source.get_step_range()
        self.step = None
        self.add_trait('step', Range(*rng))
        if self._step_editor.high == 0:
            self._step_editor.high = int(rng[1])

    def _step_changed(self, old, new):
        self._source.set_step(self.step)
        self._viewer.set_source_filename(self._source.filename)


class ViewerGUI(HasTraits):

    scene = Instance(MlabSceneModel, ())

    viewer = Instance(Viewer)
    set_step = Instance(SetStep)

##     @on_trait_change('scene.interactor')
##     def _post_init(self, name, old, new):
##         print 111
##         if self.scene.renderer:
##             pass

    def __init__(self, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0),
                 **traits):
        HasTraits.__init__(self, **traits)
        scene = self.scene.scene

        scene.foreground = fgcolor
        scene.background = bgcolor

    def show(self):
        """Show itself."""
        
        self.configure_traits()

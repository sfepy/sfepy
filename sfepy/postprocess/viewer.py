from sfepy.base.base import *
from sfepy.base.la import cycle

try:
    from enthought.mayavi import mlab
except ImportError:
    mlab = None

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
    
def add_glyphs(obj, position, bbox, rel_scaling=None,
               scale_factor='auto', clamping=False, color=None):

    glyphs = mlab.pipeline.glyph(obj, mode='2darrow', scale_mode='vector',
                                 color=color, opacity=1.0) 
    if scale_factor == 'auto':
        rng = glyphs.glyph.glyph.range
        delta = rng[1] - rng[0]
        dx = nm.max((bbox[1::2] - bbox[:-1:2]))
        if rel_scaling is None:
            rel_scaling = 0.02 # -> delta fits 50x into dx.
        scale_factor = rel_scaling * dx / delta

    glyphs.glyph.color_mode = 'color_by_vector'
    glyphs.glyph.scale_mode = 'scale_by_vector'
    glyphs.glyph.glyph.clamping = clamping
    glyphs.glyph.glyph.scale_factor = scale_factor
    glyphs.glyph.glyph_source.glyph_position = 'tail'
    glyphs.actor.actor.position = position
    return glyphs

def add_text(obj, position, text, color=(0, 0, 0)):

    t = mlab.text(x=position[0], y=position[1], text=text,
                  z=position[2], color=color, width=0.02 * len(text))
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
            
    def __call__(self, *args, **kwargs):
        """
        This is either call_mlab() or call_empty().
        """
        pass
            
    def call_empty(self, *args, **kwargs):
        pass
    
    def call_mlab(self, show=True, is_3d=False, view=None, rel_scaling=None,
                  clamping=False, layout='rowcol', fig_filename='view.png',
                  filter_names=None):
        """By default, plot all found data."""
        if filter_names is None:
            filter_names = ['node_groups', 'mat_id']

        mlab.options.offscreen = self.offscreen
        if layout == 'rowcol':
            size = (800, 600)
        elif layout == 'row':
            size = (1000, 600)
        elif layout == 'col':
            size = (600, 1000)
        else:
            size = (600, 800)
        scene = mlab.figure(bgcolor=(1,1,1), fgcolor=(0, 0, 0), size=size)

        source = mlab.pipeline.open(self.filename)
        bbox = nm.array(source.reader.unstructured_grid_output.bounds)
        dx = 1.1 * (bbox[1::2] - bbox[:-1:2])
        view3d = abs(dx[2]) > (10.0 * nm.finfo(nm.float64).eps)
        
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

        names = p_names + c_names
        names = [ii for ii in names if ii[2] not in filter_names]
        
        n_data = len(names)
        n_row, n_col = get_position_counts(n_data, layout)

        if c_names:
            ctp = mlab.pipeline.cell_to_point_data(source)

        if layout[:3] == 'col':
            iterator = enumerate(cycle((n_col, n_row)))
        else:
            iterator = enumerate(cycle((n_row, n_col)))

        for ii, (ir, ic) in iterator:
            if layout[:3] == 'col':
                ir, ic = ic, ir
            if ii == n_data: break
            family, kind, name = names[ii]

            position = nm.array([dx[0] * ic, dx[1] * (n_row - ir - 1), 0])
            position[:2] -= bbox[:2]
            
            output(family, kind, name, position)
            if kind == 'scalars':
                active = mlab.pipeline.set_active_attribute(source)
#                active.point_scalars_name = name
                setattr(active, '%s_%s_name' % (family, kind), name)

                if is_3d:
                    scp = add_scalar_cut_plane(active,
                                               position, [1, 0, 0],
                                               opacity=0.5)
                    scp = add_scalar_cut_plane(active,
                                               position, [0, 1, 0],
                                               opacity=0.5 )
                    scp = add_scalar_cut_plane(active,
                                               position, [0, 0, 1],
                                               opacity=0.5 )
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

            elif kind == 'tensors':
                if family == 'point':
                    active = mlab.pipeline.set_active_attribute(source)
                else:
                    active = mlab.pipeline.set_active_attribute(ctp)
                active.point_tensors_name = name

                active = mlab.pipeline.extract_tensor_components(active)
                if is_3d:
                    scp = add_scalar_cut_plane(active,
                                               position, [1, 0, 0],
                                               opacity=0.5)
                    scp = add_scalar_cut_plane(active,
                                               position, [0, 1, 0],
                                               opacity=0.5 )
                    scp = add_scalar_cut_plane(active,
                                               position, [0, 0, 1],
                                               opacity=0.5 )
                else:
                    surf = add_surf(active, position)

            else:
                raise ValueError('bad kind! (%s)' % kind)

            position[2] = 0.5 * dx[2]
            text = add_text(active, position, name)

            scene.scene.reset_zoom()

        scene.scene.reset_zoom()
        scene.scene.camera.zoom(1.0)

        if view is None:
            if is_3d or view3d:
                mlab.view(45, 45)
            else:
                mlab.view(0, 0)
        else:
            view = tuple(float(ii) for ii in view.split(','))
            mlab.view(*view)

        if self.auto_screenshot:
            name = os.path.join(self.output_dir, fig_filename)
            output('saving %s...' % name)
            scene.scene.save(name)
            output('...done')

        if show:
            mlab.show()

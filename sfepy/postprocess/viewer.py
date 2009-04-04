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
    
def add_glyphs(obj, position, scale_factor='auto', color=None):
    glyphs = mlab.pipeline.glyph(obj, mode='2darrow', scale_mode='vector',
                                 color=color, opacity=1.0) 
    if scale_factor == 'auto':
        rng = glyphs.glyph.glyph.range
        delta = rng[1] - rng[0]
        scale_factor = 0.01 * delta

    glyphs.glyph.color_mode = 'color_by_vector'
    glyphs.glyph.scale_mode = 'scale_by_vector'
#    glyphs.glyph.glyph.clamping = False
    glyphs.glyph.glyph.scale_factor = scale_factor
    glyphs.glyph.glyph_source.glyph_position = 'tail'
    glyphs.actor.actor.position = position
    return glyphs

def add_text(obj, position, text, color=(0, 0, 0)):

    t = mlab.text(x=position[0], y=position[1], text=text,
                  z=position[2], color=color, width=0.02 * len(text))
    return t

class Viewer(Struct):
    def __init__(self, filename, output_dir='.', offscreen=False,
                 auto_screenshot=True):
        Struct.__init__(self,
                        filename = filename,
                        output_dir = output_dir,
                        offscreen = offscreen,
                        auto_screenshot = auto_screenshot,
                        mlab = mlab)

    def __call__(self, show=True, is_3d=False):
        """By default, plot all found data."""
        mlab.options.offscreen = self.offscreen
        scene = mlab.figure(bgcolor=(1,1,1), size=(600,800))

        source = mlab.pipeline.open(self.filename)

        scalar_names = sorted( source._point_scalars_list[:-1] )
        vector_names = sorted( source._point_vectors_list[:-1] )
        names = [['scalar', name] for name in scalar_names]
        names += [['vector', name] for name in vector_names]
        
        output(scalar_names)
        output(vector_names)

        n_scalar = len(scalar_names)
        n_vector = len(vector_names)

        n_data = n_scalar + n_vector
        n_col = min(5.0, nm.fix(nm.sqrt(n_data)))
        n_row = int(nm.ceil(n_data / n_col))
        n_col = int(n_col)

        for ii, (ir, ic) in enumerate(cycle((n_row, n_col))):
            if ii == n_data: break
            kind, name = names[ii]
            
            position = [ic, n_row - ir - 1, 0]
            print position
            field_chooser = mlab.pipeline.set_active_attribute(source)
            if kind == 'scalar':
                field_chooser.point_scalars_name = name

                if is_3d:
                    scp = add_scalar_cut_plane(field_chooser,
                                               position, [1, 0, 0],
                                               opacity=0.5)
                    scp = add_scalar_cut_plane(field_chooser,
                                               position, [0, 1, 0],
                                               opacity=0.5 )
                    scp = add_scalar_cut_plane(field_chooser,
                                               position, [0, 0, 1],
                                               opacity=0.5 )
                else:
                    surf = add_surf(field_chooser, position)
                
            elif kind == 'vector':
                field_chooser.point_vectors_name = name

                glyphs = add_glyphs(field_chooser, position)

            else:
                raise ValueError('bad kind! (%s)' % kind)

            scene.scene.reset_zoom()

        scene.scene.reset_zoom()
        scene.scene.camera.zoom(1.0)

        if is_3d:
            mlab.view(45, 45)
        else:
            mlab.view(0, 0)

        if self.auto_screenshot:
            scene.scene.save(os.path.join(self.output_dir, 'view.png'))

        if show:
            mlab.show()

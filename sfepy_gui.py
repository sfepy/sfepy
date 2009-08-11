#!/usr/bin/env python
import sys
import os
import numpy as np

from enthought.traits.api import HasTraits, Range, Instance, File, Str,\
        Directory, Button, on_trait_change
from enthought.traits.ui.api import View, Item, Group, TextEditor, \
     FileEditor, DirectoryEditor
from enthought.tvtk.pyface.scene_editor import SceneEditor
from enthought.mayavi.tools.mlab_scene_model import MlabSceneModel
from enthought.mayavi.core.pipeline_base import PipelineBase
from enthought.mayavi.core.ui.mayavi_scene import MayaviScene

from sfepy.applications import pde_solve
from sfepy.fem import ProblemDefinition
from sfepy.postprocess import Viewer

class SfePyGUI(HasTraits):
    """Mayavi2-based GUI."""
    scene = Instance(MlabSceneModel, ())
    run = Button('run')

    input_filename = File('', label='input file name')
    output_dir = Directory(os.path.join(os.getcwd(), 'output'),
                           label='output directory')
    problem = Instance(ProblemDefinition)

    @on_trait_change('input_filename')
    def solve(self):
#        import pdb; pdb.set_trace()
        if not self.scene._renderer: return
        if not os.path.isfile(self.input_filename): return

        aux = pde_solve(self.input_filename, output_dir=self.output_dir)
        self.problem, self.vec, self.data = aux
        print self.vec
        self.problem, self.vec, self.data = aux
        
    def _run_fired(self):
        self.solve()

    @on_trait_change('problem')
    def view_results(self):
        self.scene.mlab.get_engine().scenes[0].scene.foreground = (0, 0, 0)
        self.scene.mlab.get_engine().scenes[0].scene.background = (1, 1, 1)
        self.scene.mlab.clf()

        view = Viewer(self.problem.get_output_name())
        view(scene=self.scene)

    # The layout of the dialog created.
    view = View(
        Item('scene', editor=SceneEditor(scene_class=MayaviScene), 
             height=600, width=800,
             show_label=False), 
        Group('_',
              Item('input_filename'),
              Item('output_dir'),
              '_'),
        Item('run', springy=True, height=10, width=20, show_label=False),
        resizable=True,
    )

gui = SfePyGUI()
if len(sys.argv) > 1:
    gui.input_filename = sys.argv[1]
gui.configure_traits()

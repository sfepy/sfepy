#!/usr/bin/env python
import sys
import os
from threading import Thread
from time import sleep
import numpy as np

from enthought.traits.api \
     import HasTraits, Instance, Int, File, Directory, Button, Code, \
     on_trait_change
from enthought.traits.ui.api \
     import View, Item, Group, CodeEditor, VSplit, HGroup, spring
from enthought.tvtk.pyface.scene_editor import SceneEditor
from enthought.pyface.timer.api import do_later
from enthought.mayavi.tools.mlab_scene_model import MlabSceneModel
from enthought.mayavi.core.ui.mayavi_scene import MayaviScene

from sfepy.applications import pde_solve
from sfepy.fem import ProblemDefinition
from sfepy.postprocess import Viewer

def assign_solution_to_gui(gui, sol):
    gui.problem, gui.vec, gui.data = sol

class SolveThread(Thread):
    def run(self):
        sol = pde_solve(self.gui.input_filename,
                        output_dir=self.gui.output_dir)
        do_later(assign_solution_to_gui, self.gui, sol)

class SfePyGUI(HasTraits):
    """Mayavi2-based GUI."""
    scene = Instance(MlabSceneModel, ())
    button_run = Button('run')
    button_clear = Button('clear output')

    input_filename = File('', label='input file name')
    output_dir = Directory(os.path.join(os.getcwd(), 'output'),
                           label='output directory')

    log = Code('')
    selected_line = Int(0)
    
    problem = Instance(ProblemDefinition)
    vec = Instance(np.ndarray)
    data = Instance(dict)
    solve_thread = Instance(SolveThread)

    @on_trait_change('input_filename')
    def solve(self):
        if not self.scene._renderer: return
        if not os.path.isfile(self.input_filename): return

        if not (self.solve_thread and self.solve_thread.isAlive()):
            self.solve_thread = SolveThread()
            self.solve_thread.gui = self
            self.solve_thread.start()
        
    def _button_run_fired(self):
        self.solve()

    def _button_clear_fired(self):
        self.log = ''

    @on_trait_change('problem')
    def view_results(self):
        self.scene.mlab.get_engine().scenes[0].scene.foreground = (0, 0, 0)
        self.scene.mlab.get_engine().scenes[0].scene.background = (1, 1, 1)
        self.scene.mlab.clf()

        view = Viewer(self.problem.get_output_name())
        view(scene=self.scene)

##     def _log_changed(self, new):
##         def aux():
##             self.selected_line = new.count('\n')  # Maybe -1 or +1, not sure.
##         do_later(aux)

    def flush():
        pass

    def write(self, text):
        def aux():
            self.log += text
            self.selected_line = self.log.count('\n')
        do_later(aux)

    # The layout of the dialog created.
    view = View(
        VSplit(
            Item('scene', editor=SceneEditor(scene_class=MayaviScene), 
                 height=600, width=800, show_label=False), 
            Group(
                Item('input_filename'),
                Item('output_dir'),
            ),
            HGroup(
                spring,
                Item('button_run', show_label=False),
                Item('button_clear', show_label=False),
            ),
            Item('log', editor=CodeEditor(auto_scroll=True,
                                          selected_line='selected_line'),
                 style='readonly', width=1.0, height=100, show_label=False),
        ),
        resizable=True,
    )

if __name__ == '__main__':

    gui = SfePyGUI()
    sys.stdout = gui
    
    if len(sys.argv) > 1:
        gui.input_filename = sys.argv[1]
        do_later(gui.solve)

    gui.configure_traits()

#!/usr/bin/env python
import sys
import os
import glob
import time
from threading import Thread

from enthought.traits.api \
     import HasTraits, Instance, Int, File, Directory, Button, Code, \
     on_trait_change
from enthought.traits.ui.api \
     import View, Item, Group, CodeEditor, VSplit, HGroup, spring
from enthought.pyface.timer.api import do_later

from sfepy.applications import pde_solve
from sfepy.fem import ProblemDefinition, State
from sfepy.postprocess import Viewer, ViewerGUI

def assign_solution_to_gui(gui, sol):
    gui.problem, gui.vec = sol

class SolveThread(Thread):
    def run(self):
        sol = pde_solve(self.gui.input_filename,
                        output_dir=self.gui.output_dir)
        do_later(assign_solution_to_gui, self.gui, sol)

class SfePyGUI(HasTraits):
    """Mayavi2-based GUI."""
    viewer_gui = Instance(ViewerGUI)
    
    button_run = Button('run')
    button_clear = Button('clear output')

    input_filename = File('', label='input file name')
    output_dir = Directory(os.path.join(os.getcwd(), 'output'),
                           label='output directory')

    log = Code('')
    selected_line = Int(0)
    
    problem = Instance(ProblemDefinition)
    vec = Instance(State)
    solve_thread = Instance(SolveThread)

    time_tag = 0.0
    text_buffer = ''

    @on_trait_change('input_filename')
    def solve(self):
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
        ts = self.problem.ts
        if ts.n_step == 1:
            output_name = self.problem.get_output_name()
        else:
            aux = self.problem.get_output_name(suffix=ts.n_digit * '?')
            output_name = sorted(glob.glob(aux))

        view = Viewer(output_name)
        self.viewer_gui = view()

    def _log_changed(self, new):
        self.selected_line = new.count('\n')  # Maybe -1 or +1, not sure.

    def flush():
        pass

    def write(self, text):
        def aux():
            self.log += self.text_buffer
            self.text_buffer = ''

        self.text_buffer += text

        tt = time.clock()
        if tt > (self.time_tag + 2.0):
            do_later(aux)
            self.time_tag = tt

    # The layout of the dialog created.
    view = View(
        VSplit(
            Item('viewer_gui', defined_when='viewer_gui is not None',
                 height=600, width=800, show_label=False,
                 style='custom'),
            Group(
                Item('input_filename'),
                Item('output_dir'),
            ),
            HGroup(
                spring,
                Item('button_run', show_label=False),
                Item('button_clear', show_label=False),
            ),
            Item('log', editor=CodeEditor(auto_scroll=False,
                                          selected_line='selected_line'),
                 style='readonly', width=1.0, height=600, show_label=False),
        ),
        resizable=True,
        width=800,
    )

if __name__ == '__main__':

    gui = SfePyGUI()
    sys.stdout = gui
    
    if len(sys.argv) > 1:
        gui.input_filename = sys.argv[1]

    gui.configure_traits()

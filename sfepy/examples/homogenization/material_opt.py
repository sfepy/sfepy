#!/usr/bin/env python
"""
See the :ref:`sec-mat_optim` tutorial for a comprehensive description of this
example.
"""
import sys
sys.path.append('.')

import numpy as nm
from scipy.optimize import fmin_tnc

import sfepy
from sfepy.base.base import Struct
from sfepy.base.log import Log

class MaterialOptimizer:

    @staticmethod
    def create_app(filename, is_homog=False, **kwargs):
        from sfepy.base.conf import ProblemConf, get_standard_keywords
        from sfepy.homogenization.homogen_app import HomogenizationApp
        from sfepy.applications import PDESolverApp

        required, other = get_standard_keywords()
        if is_homog:
            required.remove('equations')

        conf = ProblemConf.from_file(filename, required, other,
                                     define_args=kwargs)
        options = Struct(output_filename_trunk=None,
                         save_ebc=False,
                         save_ebc_nodes=False,
                         save_regions=False,
                         save_regions_as_groups=False,
                         solve_not=False)

        if is_homog:
            app = HomogenizationApp(conf, options, 'material_opt_micro:')

        else:
            app = PDESolverApp(conf, options, 'material_opt_macro:')

        app.conf.opt_data = {}
        opts = conf.options
        if hasattr(opts, 'parametric_hook'):  # Parametric study.
            parametric_hook = conf.get_function(opts.parametric_hook)
            app.parametrize(parametric_hook)

        return app

    def x_norm2real(self, x):
        return x * (self.x_U - self.x_L) + self.x_L

    def x_real2norm(self, x):
        return (x - self.x_L) / (self.x_U - self.x_L)

    def __init__(self, macro_fn, micro_fn, x0, x_L, x_U, exp_data):
        self.macro_app = self.create_app(macro_fn, is_homog=False, is_opt=True)
        self.micro_app = self.create_app(micro_fn, is_homog=True, is_opt=True)
        self.x_L = nm.array(x_L)
        self.x_U = nm.array(x_U)
        self.x0 = self.x_real2norm(nm.array(x0))
        self.x = []
        self.eval_f = []
        self.exp_data = exp_data

    @staticmethod
    def rotate_mat(D, angle):
        s = nm.sin(angle)
        c = nm.cos(angle)
        s2 = s**2
        c2 = c**2
        sc = s * c
        T = nm.array([[c2, 0, s2, 0, 2*sc,0],
                      [0, 1, 0, 0, 0, 0],
                      [s2, 0, c2, 0, -2*sc, 0],
                      [0, 0, 0, c, 0, -s],
                      [-sc, 0, sc, 0, c2 - s2, 0],
                      [0, 0, 0, s, 0, c]])

        return nm.dot(nm.dot(T, D), T.T)

    def matopt_eval(self, x):
        mic_od = self.micro_app.conf.opt_data
        mac_od = self.macro_app.conf.opt_data

        mic_od['coefs'] = {}
        mic_od['mat_params'] = x
        self.micro_app()

        D = mic_od['D_homog']
        val = 0.0
        aux = []
        for phi, exp_k in self.exp_data:
            print('phi = %d' % phi)

            mac_od['D_homog'] = self.rotate_mat(D, nm.deg2rad(phi))
            self.macro_app()

            comp_k = mac_od['k']
            val += (1.0 - comp_k / exp_k)**2
            aux.append((comp_k, exp_k))

        val = nm.sqrt(val)
        self.x.append(x)
        self.eval_f.append(val)

        return val

    def iter_step(self, x, first_step=False):
        if first_step:
            self.log = Log([['O'], ['E_f', 'E_m'], ['v_f', 'v_m']],
                ylabels=['Obj. fun.', "Young's modulus", "Poisson's ratio"],
                xlabels=['iter', 'iter', 'iter'],
                aggregate=0)
            self.istep = 0
            self.log(0.5, x[0], x[2], x[1], x[3],
                x=[0, 0, 0, 0])
        else:
            self.istep += 1
            self.log(self.eval_f[-1], x[0], x[2], x[1], x[3],
                x=(self.istep,)*4)

    def material_optimize(self):
        x0 = self.x0
        bnds = zip(self.x_real2norm(self.x_L), self.x_real2norm(self.x_U))
        feval = lambda x: self.matopt_eval(self.x_norm2real(x))
        istep = lambda x: self.iter_step(self.x_norm2real(x))
        self.iter_step(self.x_norm2real(x0), first_step=True)

        print('>>> material optimization START <<<')
        xopt = fmin_tnc(feval, x0, approx_grad=True, bounds=list(bnds),
                        xtol=1e-3, callback=istep)
        print('>>> material optimization FINISHED <<<')

        self.log(finished=True)
        return self.x_norm2real(xopt[0])

def main():
    srcdir = sfepy.base_dir + '/examples/homogenization/'
    micro_filename = srcdir + 'homogenization_opt.py'
    macro_filename = srcdir + 'linear_elasticity_opt.py'

    exp_data = zip([0, 30, 60, 90], [1051140., 197330., 101226., 95474.])
    mo = MaterialOptimizer(macro_filename, micro_filename,
                           [160.e9, 0.25, 5.e9, 0.45],
                           [120e9, 0.2, 2e9, 0.2],
                           [200e9, 0.45, 8e9, 0.45],
                           list(exp_data))

    optim_par = mo.material_optimize()
    print('optimized parameters: ', optim_par)

if __name__ == '__main__':
    main()

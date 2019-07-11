r"""
DG FEM convergence test for 1D diffusion only problem
"""
import numpy as nm
import sympy as sm


from examples.dg.example_dg_common import get_1Dmesh_hook
from examples.dg.conv_tests.fem_conv_test import SimpleExpression, define_functions, eval_expr
from os.path import join as pjoin

from my_utils.visualizer import reconstruct_legendre_dofs


def define(filename_mesh=None, approx_order=1, Cw=100, diffusion_coef=0.0001, use_symbolic=False):

    if filename_mesh is None:
        filename_mesh = get_1Dmesh_hook(0, 1, 2)


    materials = {
        'D': ({'val': [diffusion_coef], '.Cw': Cw},),
        'g': 'source_fun'
    }

    regions = {
        'Omega' : 'all',
        'Gamma' : ('vertices of surface', 'facet'),
        'left': ('vertices in x == 0', 'vertex'),
        'right': ('vertices in x == 1', 'vertex')
    }

    fields = {
        'f': ('real', 'scalar', 'Omega', str(approx_order) + 'd', 'DG', 'legendre')
    }

    variables = {
        'u' : ('unknown field', 'f', 0),
        'v' : ('test field',    'f', 'u'),
    }

    # dgebcs = {
    #     'u_left': ('left', {'u.all': "bc_fun", 'grad.u.all': "bc_fun"}),
    #     'u_right': ('right', {'u.all': "bc_fun", 'grad.u.all': "bc_fun"}),
    # }

    dgepbc_1 = {
        'name'  : 'u_rl',
        'region': ['right', 'left'],
        'dofs': {'u.all': 'u.all'},
        'match': 'match_y_line',
    }

    integrals = {
        'i' : 2 * approx_order,
    }

    equations = {
        'Temperature': " - dw_laplace.i.Omega(D.val, v, u) " +
                       " + dw_dg_diffusion_flux.i.Omega(D.val, u, v)" +
                       " + dw_dg_diffusion_flux.i.Omega(D.val, v, u)" +
                       " - " + str(diffusion_coef) + "* dw_dg_interior_penal.i.Omega(D.Cw, v, u)" +
                       " + dw_volume_lvf.i.Omega(g.val, v) = 0"
    }

    solvers = {
        'ls' : ('ls.auto_direct', {}),
        'newton' : ('nls.newton',
                    {'i_max'      : 1,
                     'eps_a'      : 1e-10,
        }),
    }

    options = {
        'nls' : 'newton',
        'ls' : 'ls',
        'output_format'   : 'msh',

    }

    functions = {}
    def local_register_function(fun):
        try:
            functions.update({fun.__name__: (fun,)})

        except AttributeError:  # Already a sfepy Function.
            fun = fun.function
            functions.update({fun.__name__: (fun,)})

        return fun

    if use_symbolic:
        estr = '-(y ** 2 - y) * sin(2 * pi * x)'
        expression = SimpleExpression(estr)
        expr = expression.define_expression()
        bc_fun, source_fun, sol_fun = expression.define_functions(expr)

        bc_fun = local_register_function(bc_fun)
        source_fun = local_register_function(source_fun)
        sol_fun = local_register_function(sol_fun)

    else:
        @local_register_function
        def bc_fun(ts, coors, bc, problem):
            x = coors[..., 0]
            res = nm.zeros(nm.shape(x))
            pi = nm.pi

            if bc.diff == 1:
                if "left" in bc.name:
                    res[:] = 2*pi
                elif "right" in bc.name:
                    res[:] = 2*pi
            else:
                # bc are zero
                pass
            return res

        @local_register_function
        def source_fun(ts, coors, mode="qp", **kwargs):
            # t = ts.dt * ts.step
            eps = diffusion_coef
            sin = nm.sin
            cos = nm.cos
            pi = nm.pi
            if mode == "qp":
                x = coors[..., 0]
                res = 4*pi**2*eps*sin(2*pi*x)
                return {"val": res[..., None, None]}

        def analytic_sol(coors, t=0):
            x = coors[..., 0]
            sin = nm.sin
            pi = nm.pi
            res = sin(2*pi*x)
            return res

        @local_register_function
        def sol_fun(ts, coors, mode="qp", **kwargs):
            t = ts.time
            if mode == "qp":
                return {"u": analytic_sol(coors, t)[..., None, None]}

    return locals()

def main():
    import sys
    import time
    from sfepy.base.ioutils import ensure_path
    from sfepy.base.base import output
    from sfepy.base.conf import ProblemConf
    from sfepy.discrete import (Integral, Integrals, Material, Problem)
    from sfepy.discrete.common.mappings import get_jacobian
    from sfepy.discrete.dg.my_utils.plot_1D_dg import clear_folder
    from matplotlib import pyplot as plt

    base_output_folder = "output\\diff1D\\"

    n_nod = 5
    n_refine = 6
    orders = [1, 2, 3, 4, 5]
    #orders = [1]

    mod = sys.modules[__name__]

    sol_fig, axs = plt.subplots(len(orders), n_refine, figsize=(18, 10))

    results = []
    for ir, refine in enumerate(range(n_refine)):
        gen_mesh = get_1Dmesh_hook(0, 1, n_nod)
        for io, order in enumerate(orders):
            output('n_nod:', n_nod, 'order:', order)

            tt = time.clock()
            conf = ProblemConf.from_dict(define(gen_mesh, order), mod)
            pb = Problem.from_conf(conf)
            sol = pb.solve()
            elapsed = time.clock() - tt

            output_folder = pjoin(base_output_folder, "h" + str(n_nod - 1))
            output_folder = pjoin(output_folder, "o" + str(order))

            output_format = pjoin(output_folder, "sol-h{:02d}o{:02d}.*.{}".format(n_nod, order, "vtk"))
            output("Output set to {}, clearing ...".format(output_format))

            clear_folder(output_format, confirm=False)
            ensure_path(output_format)

            pb.save_state(output_format.replace("*", "0"), state=sol)

            idiff = Integral('idiff', 20)

            qps = pb.fields["f"].mapping.get_physical_qps(idiff.get_qp("1_2")[0])
            fqps = qps.flatten()

            num_qp = pb.evaluate(
                'ev_volume_integrate.idiff.Omega(u)',
                integrals=Integrals([idiff]), mode='qp', copy_materials=False,
            )

            aux = Material('aux', function=conf.sol_fun)
            ana_qp = pb.evaluate(
                'ev_volume_integrate_mat.idiff.Omega(aux.u, u)',
                aux=aux, integrals=Integrals([idiff]), mode='qp',
                copy_materials=False,
            )

            field = pb.fields['f']
            det = get_jacobian(field, idiff)

            diff_l2 = nm.sqrt((((num_qp - ana_qp)**2) * det).sum())
            ana_l2 = nm.sqrt((((ana_qp)**2) * det).sum())
            error = diff_l2 / ana_l2

            n_dof = field.n_nod
            sol_fig.suptitle(
                "Numerical and exact solutions, Cw: {}, diffusion: {}".format(conf.Cw, conf.diffusion_coef))

            # from sfepy.discrete.dg.my_utils.visualizer import plot_1D_legendre_dofs
            coors = pb.domain.mesh.coors
            u = pb.fields["f"].unravel_sol(sol.vec)
            uu, xx = reconstruct_legendre_dofs(coors, None, u.swapaxes(0, 1)[:, :, None])

            ax = axs[order-1][refine]

            xs = nm.linspace(nm.min(0), nm.max(1), 500)[:, None]
            ax.set_title("o: {}, h: {}".format(order, n_nod - 1))
            ax.plot(xs, conf.analytic_sol(xs), label="fun-ex", color="grey")
            ax.plot(xx[:, 0], uu[:, 0, 0], alpha=.5)
            ax.plot(fqps, ana_qp.flatten(), "--", color="grey")
            ax.plot(fqps, num_qp.flatten())
            ax2 = ax.twinx()
            ax2.plot(fqps, nm.abs(num_qp.flatten() - ana_qp.flatten()), color="red")
            if io < len(orders) - 1:
                ax.set_xticks([])
            # if ir > 0:
            #     ax.set_yticks([])

            result = (n_nod-1, order, n_dof, ana_l2, diff_l2, error, elapsed)
            results.append(result)

        n_nod = n_nod + n_nod - 1

    sol_fig.savefig("per-err-sol-i20cw{}_d{}.tif".format(conf.Cw, conf.diffusion_coef), dpi=100)
    results = nm.array(results)
    output(results)

    conv_fig = plt.figure()
    conv_fig.suptitle("Convergences by order, Cw: {}, diffusion: {}".format(conf.Cw, conf.diffusion_coef))
    for o in orders:
        curr_results = results[results[:, 1] == o]
        co = plt.loglog(1/curr_results[:, 0], curr_results[:, 4], 'o', label=str(o))[0].get_color()
        plt.loglog(1/curr_results[:, 0], curr_results[:, 4], color=co)
        plt.grid()
        plt.xlabel("h")
        plt.ylabel("L^2 error")
    plt.legend()
    conv_fig.savefig("per-conv-i20cw{}_d{}.tif".format(conf.Cw, conf.diffusion_coef), dpi=200)


    plt.show()

if __name__ == '__main__':
    main()

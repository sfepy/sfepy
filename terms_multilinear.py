import opt_einsum as oe

from sfepy.base.base import Struct
from sfepy.terms.terms import Term
from sfepy.terms import register_term

class ETermBase(Struct):

    def einsum(self, sexpr, *args, diff_var=None):
        dc_type = self.get_dof_conn_type()

        vvar = args[0]
        uvars = args[1:]
        vg, _ = self.get_mapping(vvar)

        dets = vg.det
        qsb = vg.bf
        qsbg = vg.bfg

        expr = ['cqab']
        eargs = [dets]

        eins = sexpr.split(',')
        letters = 'defgh'
        out_expr = 'c'
        for ii, ein in enumerate(eins):
            arg = args[ii]
            if arg.is_virtual():
                if arg.n_components == 1:
                    if '.' in ein: # derivative
                        expr.append('cq{}{}'.format(ein[2], letters[ii]))
                        eargs.append(qsbg)
                        out_expr += letters[ii]

                    else:
                        raise NotImplementedError
                else:
                    raise NotImplementedError

            else:
                if arg.n_components == 1:
                    if '.' in ein: # derivative
                        expr.append('cq{}{}'.format(ein[2], letters[ii]))
                        eargs.append(qsbg)
                    else:
                        raise NotImplementedError

                    if diff_var != arg.name:
                        # Assumes no E(P)BCs are present!
                        adc = arg.get_dof_conn(dc_type)
                        dofs = arg()[adc]

                        expr.append('c{}'.format(letters[ii]))
                        eargs.append(dofs)

                    else:
                        out_expr += letters[ii]

                else:
                    raise NotImplementedError

        self.parsed_expression = ','.join(expr) + '->' + out_expr
        self.path, self.path_info = oe.contract_path(
            self.parsed_expression, *eargs
        )
        if diff_var is not None:
            def eval_einsum(out):
                return oe.contract(self.parsed_expression, *eargs,
                                   out=out[:, 0, ...], optimize=self.path)

        else:
            def eval_einsum(out):
                return oe.contract(self.parsed_expression, *eargs,
                                   out=out[:, 0, :, 0], optimize=self.path)

        return eval_einsum

    @staticmethod
    def function(out, eval_einsum):
        eval_einsum(out)
        return 0

    def get_fargs(self, *args, **kwargs):
        eval_einsum = self.expression(*args, **kwargs)

        return eval_einsum,

class ELaplaceTerm(ETermBase, Term):
    name = 'dw_elaplace'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'parameter_1', 'parameter_2'))
    arg_shapes = [{'opt_material' : '1, 1', 'virtual' : (1, 'state'),
                   'state' : 1, 'parameter_1' : 1, 'parameter_2' : 1},
                  {'opt_material' : None}]
    modes = ('weak', 'eval')

    def expression(self, mat, virtual, state, mode=None, term_mode=None,
                   diff_var=None, **kwargs):
        """
        diff_var not needed here(?), but Term passes it in *args.
        """
        if mat is None:
            expr = self.einsum('0.j,0.j', virtual, state,
                               diff_var=diff_var)

        else:
            expr = self.einsum('jk,0.j,0.k', mat, virtual, state,
                               diff_var=diff_var)

        return expr

register_term(ELaplaceTerm)

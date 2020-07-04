import numpy as nm
import opt_einsum as oe

from sfepy.base.base import Struct
from sfepy.terms.terms import Term
from sfepy.terms import register_term

class ETermBase(Struct):

    def einsum(self, sexpr, *args, diff_var=None):
        dc_type = self.get_dof_conn_type()

        vvar = args[0]
        uvars = args[1:]

        if diff_var is not None:
            n_add = len([var.name for var in uvars if var.name == diff_var])

        else:
            n_add = 1

        vg, _ = self.get_mapping(vvar)

        dets = vg.det
        qsb = vg.bf
        qsbg = vg.bfg

        exprs = [['cqab']] * n_add
        eargss = [[dets]] * n_add
        def append_all(seqs, item):
            for seq in seqs:
                seq.append(item)

        eins = sexpr.split(',')
        letters = 'defgh'
        out_expr = ['c'] * n_add
        ia = 0
        for ii, ein in enumerate(eins):
            arg = args[ii]

            if arg.is_virtual():
                if arg.n_components == 1:
                    if '.' in ein: # derivative
                        append_all(exprs, 'cq{}{}'.format(ein[2], letters[ii]))
                        append_all(eargss, qsbg)
                        for iia in range(n_add):
                            out_expr[iia] += letters[ii]

                    else:
                        raise NotImplementedError
                else:
                    raise NotImplementedError

            else:
                if arg.n_components == 1:
                    if '.' in ein: # derivative
                        eterm = 'cq{}{}'.format(ein[2], letters[ii])
                        earg = qsbg
                    else:
                        raise NotImplementedError

                    append_all(exprs, eterm)
                    append_all(eargss, earg)
                    if (diff_var != arg.name) or (n_add > 1):
                        # Assumes no E(P)BCs are present!
                        adc = arg.get_dof_conn(dc_type)
                        dofs = arg()[adc]

                        determ = 'c{}'.format(letters[ii])

                    if (diff_var != arg.name):
                        append_all(exprs, determ)
                        append_all(eargss, dofs)

                    else:
                        for iia in range(n_add):
                            if iia != ia:
                                exprs[iia].append(determ)

                        out_expr[ia] += letters[ii]
                        ia += 1

                else:
                    raise NotImplementedError

        self.parsed_expressions = [','.join(exprs[ia]) + '->' + out_expr[ia]
                                   for ia in range(n_add)]
        self.paths, self.path_infos = zip(*[oe.contract_path(
            self.parsed_expressions[ia], *eargss[ia], optimize='auto',
        ) for ia in range(n_add)])
        print(self.parsed_expressions)

        if n_add == 1:
            if diff_var is not None:
                def eval_einsum(out):
                    oe.contract(self.parsed_expressions[0], *eargss[0],
                                out=out[:, 0, ...],
                                optimize=self.paths[0])

            else:
                def eval_einsum(out):
                    oe.contract(self.parsed_expressions[0], *eargss[0],
                                out=out[:, 0, :, 0],
                                optimize=self.paths[0])

        else:
            if diff_var is not None:
                def eval_einsum(out):
                    aux = nm.empty_like(out)
                    for ia in range(n_add):
                        oe.contract(self.parsed_expressions[ia], *eargss[ia],
                                    out=aux[:, 0, ...],
                                    optimize=self.paths[ia])
                        out[:] += aux

            else:
                def eval_einsum(out):
                    aux = nm.empty_like(out)
                    for ia in range(n_add):
                        oe.contract(self.parsed_expressions[ia], *eargss[ia],
                                    out=out[:, 0, :, 0],
                                    optimize=self.paths[ia])
                        out[:] += aux

        return eval_einsum

    @staticmethod
    def function(out, eval_einsum):
        eval_einsum(out)
        return 0

    def get_fargs(self, *args, **kwargs):
        # This should be called on construction?
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

class EConvectTerm(ETermBase, Term):
    name = 'dw_econvect'
    arg_types = ('virtual', 'state')
    arg_shapes = {'virtual' : ('D', 'state'), 'state' : 'D'}

    def expression(self, mat, virtual, state, mode=None, term_mode=None,
                   diff_var=None, **kwargs):
        expr = self.einsum('i,j.i,j', virtual, state, state,
                           diff_var=diff_var)

        return expr

register_term(EConvectTerm)

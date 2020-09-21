import numpy as nm
import opt_einsum as oe

from sfepy.base.base import Struct
from sfepy.base.timing import Timer
from sfepy.terms.terms import Term
from sfepy.terms import register_term

def _get_char_map(c1, c2):
    mm = {}
    for ic, char in enumerate(c1):
        if char in mm:
            print(char, '->eq?', mm[char], c2[ic])
            if mm[char] != c2[ic]:
                mm[char] += c2[ic]
        else:
            mm[char] = c2[ic]

    return mm

class ETermBase(Struct):

    """
    Reserved letters:

    c .. cells
    q .. quadrature points
    d-h .. DOFs axes
    r-z .. auxiliary axes
    """

    def einsum(self, sexpr, *args, diff_var=None):

        timer = Timer('')
        timer.start()

        dc_type = self.get_dof_conn_type()

        vvar = args[0]
        uvars = args[1:]

        n_elr, n_qpr, dim, n_enr, n_cr = self.get_data_shape(vvar)

        if diff_var is not None:
            n_add = len([var.name for var in uvars if var.name == diff_var])

            varc = self.get_variables(as_list=False)[diff_var]
            n_elc, n_qpc, dim, n_enc, n_cc = self.get_data_shape(varc)
            eshape = tuple([n_elr]
                           + ([n_cr] if n_cr > 1 else [])
                           + [n_enr]
                           + ([n_cc] if n_cc > 1 else [])
                           + [n_enc])

        else:
            n_add = 1
            eshape = (n_elr, n_cr, n_enr) if n_cr > 1 else (n_elr, n_enr)

        vg, _ = self.get_mapping(vvar)

        dets = vg.det
        qsb = vg.bf
        qsbg = vg.bfg

        exprs = [['cq'] for ii in range(n_add)]
        eargss = [[dets[..., 0, 0]] for ii in range(n_add)]
        def append_all(seqs, item):
            for seq in seqs:
                seq.append(item)
        used_dofs = {}

        eins = sexpr.split(',')
        letters = 'defgh'
        aux_letters = iter('rstuvwxyz')
        out_expr = ['c'] * n_add
        ia = 0

        n_cell, n_qp, dim, n_ed = qsbg.shape
        ee = nm.eye(dim)

        for ii, ein in enumerate(eins):
            arg = args[ii]

            if arg.is_virtual():
                if arg.n_components == 1:
                    if '.' in ein: # derivative
                        append_all(exprs, 'cq{}{}'.format(ein[2], letters[ii]))
                        append_all(eargss, qsbg)

                    else:
                        append_all(exprs, 'q{}{}'.format(ein[0], letters[ii]))
                        append_all(eargss, qsb[0])

                    for iia in range(n_add):
                        out_expr[iia] += letters[ii]

                else:
                    if '.' in ein: # derivative
                        raise NotImplementedError

                    else:
                        aux = next(aux_letters)
                        iin = letters[ii] # node (qs basis index)
                        iic = next(aux_letters) # component
                        append_all(exprs, 'q{}{}'.format(aux, iin))
                        append_all(exprs, '{}{}'.format(ein[0], iic))
                        append_all(eargss, qsb[0])
                        append_all(eargss, ee)

                    for iia in range(n_add):
                        out_expr[iia] += (iic + iin)

            else:
                if arg.n_components == 1:
                    if '.' in ein: # derivative
                        eterm = 'cq{}{}'.format(ein[2], letters[ii])
                        earg = qsbg

                    else:
                        eterm = 'q{}{}'.format(ein[0], letters[ii])
                        earg = qsb[0]

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
                    aux = next(aux_letters)
                    iin = letters[ii] # node (qs basis index)
                    iic = next(aux_letters) # component

                    if '.' in ein: # derivative
                        eterm = 'cq{}{}'.format(ein[2], iin)
                        earg = qsbg

                    else:
                        eterm = 'q{}{}'.format(aux, iin)
                        earg = qsb[0]

                    append_all(exprs, eterm)
                    append_all(eargss, earg)
                    if (diff_var != arg.name) or (n_add > 1):
                        dofs = used_dofs.get(arg.name)
                        if dofs is None:
                            # Assumes no E(P)BCs are present!
                            adc = arg.get_dof_conn(dc_type)
                            dofs = arg()[adc]
                            dofs.shape = (dets.shape[0], -1, qsb.shape[-1])
                            used_dofs[arg.name] = dofs

                        determ = 'c{}{}'.format(ein[0], iin)

                    if (diff_var != arg.name):
                        append_all(exprs, determ)
                        append_all(eargss, dofs)

                    else:
                        for iia in range(n_add):
                            if iia != ia:
                                exprs[iia].append(determ)
                                eargss[iia].append(dofs)

                            else:
                                eeterm = '{}{}'.format(ein[0], iic)
                                exprs[iia].append(eeterm)
                                eargss[iia].append(ee)

                        out_letters = (iic + iin)
                        out_expr[ia] += out_letters
                        ia += 1

        self.parsed_expressions = [','.join(exprs[ia]) + '->' + out_expr[ia]
                                   for ia in range(n_add)]
        print(self.parsed_expressions)
        self.paths, self.path_infos = zip(*[oe.contract_path(
            self.parsed_expressions[ia], *eargss[ia], optimize='greedy',
        ) for ia in range(n_add)])
        print(self.paths)
        print(self.path_infos)

        if n_add == 1:
            if diff_var is not None:
                def eval_einsum(out):
                    vout = out.reshape(eshape)
                    oe.contract(self.parsed_expressions[0], *eargss[0],
                                out=vout,
                                optimize=self.paths[0])

            else:
                def eval_einsum(out):
                    tt = Timer('')
                    tt.start()

                    vout = out.reshape(eshape)
                    oe.contract(self.parsed_expressions[0], *eargss[0],
                                out=vout,
                                optimize=self.paths[0])
                    # Below is faster, due to repeated 'z' (of size 1)!
                    # oe.contract('cqab,qzy,jx,cqkl,cjl,qzn,ckn->cxy',
                    #             *eargss[0],
                    #             out=vout,
                    #             optimize=self.paths[0])
                    print(tt.stop())
                    # mm = _get_char_map(self.parsed_expressions[0],
                    #                    'cqab,qzy,jx,cqkl,cjl,qzn,ckn->cxy')
                    # for key, val in mm.items():
                    #     print(key, val)
                    # print(len(mm), len(set(mm.values())))
                    # from sfepy.base.base import debug; debug()

        else:
            if diff_var is not None:
                def eval_einsum(out):

                    # mm = _get_char_map(self.parsed_expressions[0],
                    #                    'cqab,qzy,jx,cqkY,jX,qzn,ckn->cxyXY')
                    # for key, val in mm.items():
                    #     print(key, val)
                    # print(len(mm), len(set(mm.values())))

                    # mm = _get_char_map(self.parsed_expressions[1],
                    #                    'cqab,qzy,jx,cqkl,cjl,qzY,kX->cxyXY')
                    # for key, val in mm.items():
                    #     print(key, val)
                    # print(len(mm), len(set(mm.values())))
                    # from sfepy.base.base import debug; debug()
                    tt = Timer('')
                    tt.start()

                    vout = out.reshape(eshape)
                    oe.contract(self.parsed_expressions[0], *eargss[0],
                                out=vout,
                                optimize=self.paths[0])
                    aux = nm.empty_like(out)
                    vaux = out.reshape(eshape)
                    for ia in range(1, n_add):
                        oe.contract(self.parsed_expressions[ia], *eargss[ia],
                                    out=vaux,
                                    optimize=self.paths[ia])
                        out[:] += aux

                    print(tt.stop())

            else: # This never happens?
                raise RuntimeError('Impossible code path!')

        print(timer.stop())
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

class EVolumeDotTerm(ETermBase, Term):
    name = 'dw_evolume_dot'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'parameter_1', 'parameter_2'))
    arg_shapes = [{'opt_material' : '1, 1', 'virtual' : (1, 'state'),
                   'state' : 1, 'parameter_1' : 1, 'parameter_2' : 1},
                  {'opt_material' : None},
                  {'opt_material' : '1, 1', 'virtual' : ('D', 'state'),
                   'state' : 'D', 'parameter_1' : 'D', 'parameter_2' : 'D'},
                  {'opt_material' : 'D, D'},
                  {'opt_material' : None}]
    modes = ('weak', 'eval')

    def expression(self, mat, virtual, state, mode=None, term_mode=None,
                   diff_var=None, **kwargs):
        if mat is None:
            expr = self.einsum('i,i', virtual, state,
                               diff_var=diff_var)

        else:
            expr = self.einsum('ij,i,j', mat, virtual, state,
                               diff_var=diff_var)

        return expr

register_term(EVolumeDotTerm)

class EConvectTerm(ETermBase, Term):
    name = 'dw_econvect'
    arg_types = ('virtual', 'state')
    arg_shapes = {'virtual' : ('D', 'state'), 'state' : 'D'}

    def expression(self, virtual, state, mode=None, term_mode=None,
                   diff_var=None, **kwargs):
        expr = self.einsum('i,i.j,j', virtual, state, state,
                           diff_var=diff_var)

        return expr

register_term(EConvectTerm)

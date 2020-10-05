import numpy as nm
import opt_einsum as oe

from sfepy.base.base import output, Struct
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

def append_all(seqs, item, ii=None):
    if ii is None:
        for seq in seqs:
            seq.append(item)

    else:
        seqs[ii].append(item)

def get_sizes(indices, operands):
    sizes = {}

    for iis, op in zip(indices, operands):
        for ii, size in zip(iis, op.shape):
            sizes[ii] = size

    return sizes

class ExpressionBuilder(Struct):
    letters = 'defgh'
    aux_letters = iter('rstuvwxyz')

    def __init__(self, n_add, dc_type, dofs_cache):
        self.n_add = n_add
        self.subscripts = [[] for ia in range(n_add)]
        self.operands = [[] for ia in range(n_add)]
        self.operand_names = [[] for ia in range(n_add)]
        self.out_subscripts = ['c' for ia in range(n_add)]
        self.ia = 0
        self.dc_type = dc_type
        self.dofs_cache = dofs_cache

    def add_constant(self, val, name):
        append_all(self.subscripts, 'cq')
        append_all(self.operands, val)
        append_all(self.operand_names, name)

    def add_bfg(self, iin, ein, qsbg, name):
        append_all(self.subscripts, 'cq{}{}'.format(ein[2], iin))
        append_all(self.operands, qsbg)
        append_all(self.operand_names, name)

    def add_bf(self, iin, ein, qsb, name):
        append_all(self.subscripts, 'q{}'.format(iin))
        append_all(self.operands, qsb[0, :, 0])
        append_all(self.operand_names, name)

    def add_eye(self, iic, ein, eye, iia=None):
        append_all(self.subscripts, '{}{}'.format(ein[0], iic), ii=iia)
        append_all(self.operands, eye, ii=iia)
        append_all(self.operand_names, 'I', ii=iia)

    def add_arg_dofs(self, iin, ein, arg, iia=None):
        dofs = self.dofs_cache.get(arg.name)
        if dofs is None:
            # Assumes no E(P)BCs are present!
            adc = arg.get_dof_conn(self.dc_type)
            dofs = arg()[adc]
            self.dofs_cache[arg.name] = dofs

        if arg.n_components > 1:
            dofs.shape = (dofs.shape[0], arg.n_components, -1)
            term = 'c{}{}'.format(ein[0], iin)

        else:
            term = 'c{}'.format(iin)

        append_all(self.subscripts, term, ii=iia)
        append_all(self.operands, dofs, ii=iia)
        append_all(self.operand_names, arg.name + '.dofs', ii=iia)

    def add_virtual_arg(self, arg, ii, ein, qsb, qsbg):
        iin = self.letters[ii] # node (qs basis index)
        if '.' in ein: # derivative
            self.add_bfg(iin, ein, qsbg, arg.name + '.bfg')

        else:
            self.add_bf(iin, ein, qsb, arg.name + '.bf')

        out_letters = iin

        if arg.n_components > 1:
            iic = next(self.aux_letters) # component
            ee = nm.eye(arg.n_components)
            self.add_eye(iic, ein, ee)

            out_letters = iic + out_letters

        for iia in range(self.n_add):
            self.out_subscripts[iia] += out_letters

    def add_state_arg(self, arg, ii, ein, qsb, qsbg, diff_var):
        iin = self.letters[ii] # node (qs basis index)
        if '.' in ein: # derivative
            self.add_bfg(iin, ein, qsbg, arg.name + '.bfg')

        else:
            self.add_bf(iin, ein, qsb, arg.name + '.bf')

        out_letters = iin

        if (diff_var != arg.name):
            self.add_arg_dofs(iin, ein, arg)

        else:
            if arg.n_components > 1:
                iic = next(self.aux_letters) # component
                ee = nm.eye(arg.n_components)

                out_letters = iic + out_letters

            for iia in range(self.n_add):
                if iia != self.ia:
                    self.add_arg_dofs(iin, ein, arg, iia)

                elif arg.n_components > 1:
                    self.add_eye(iic, ein, ee, iia)

            self.out_subscripts[self.ia] += out_letters
            self.ia += 1

    def get_expressions(self):
        expressions = [','.join(self.subscripts[ia]) +
                       '->' +
                       self.out_subscripts[ia]
                       for ia in range(self.n_add)]
        return expressions

    def print_shapes(self):
        for ia in range(self.n_add):
            sizes = get_sizes(self.subscripts[ia], self.operands[ia])
            out_shape = tuple(sizes[ii] for ii in self.out_subscripts[ia])
            output(sizes)
            output(self.out_subscripts[ia], out_shape, '=')

            for name, ii, op in zip(self.operand_names[ia],
                                    self.subscripts[ia],
                                    self.operands[ia]):
                output('  {:10}{:8}{}'.format(name, ii, op.shape))

class ETermBase(Struct):
    """
    Reserved letters:

    c .. cells
    q .. quadrature points
    d-h .. DOFs axes
    r-z .. auxiliary axes
    """
    optimize = 'dynamic-programming'

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

        dofs_cache = {}
        self.ebuilder = ExpressionBuilder(n_add, dc_type, dofs_cache)
        self.ebuilder.add_constant(dets[..., 0, 0], 'J')

        eins = sexpr.split(',')
        for ii, ein in enumerate(eins):
            arg = args[ii]

            if arg.is_virtual():
                self.ebuilder.add_virtual_arg(arg, ii, ein, qsb, qsbg)

            else:
                ag, _ = self.get_mapping(arg)
                self.ebuilder.add_state_arg(arg, ii, ein, ag.bf, ag.bfg,
                                            diff_var)

        self.parsed_expressions = self.ebuilder.get_expressions()
        output(self.parsed_expressions)
        self.ebuilder.print_shapes()
        operands = self.ebuilder.operands
        self.paths, self.path_infos = zip(*[oe.contract_path(
            self.parsed_expressions[ia], *operands[ia],
            optimize=self.optimize,
        ) for ia in range(n_add)])
        output(self.paths)
        output(self.path_infos)

        if n_add == 1:
            if diff_var is not None:
                def eval_einsum(out):
                    tt = Timer('')
                    tt.start()

                    vout = out.reshape(eshape)
                    oe.contract(self.parsed_expressions[0], *operands[0],
                                out=vout,
                                optimize=self.paths[0])

                    output('eval_einsum 1M: {} s'.format(tt.stop()))

            else:
                def eval_einsum(out):
                    tt = Timer('')
                    tt.start()

                    vout = out.reshape(eshape)
                    oe.contract(self.parsed_expressions[0], *operands[0],
                                out=vout,
                                optimize=self.paths[0])

                    output('eval_einsum 1V: {} s'.format(tt.stop()))

        else:
            if diff_var is not None:
                def eval_einsum(out):
                    tt = Timer('')
                    tt.start()

                    vout = out.reshape(eshape)
                    oe.contract(self.parsed_expressions[0], *operands[0],
                                out=vout,
                                optimize=self.paths[0])
                    aux = nm.empty_like(out)
                    vaux = aux.reshape(eshape)
                    for ia in range(1, n_add):
                        oe.contract(self.parsed_expressions[ia], *operands[ia],
                                    out=vaux,
                                    optimize=self.paths[ia])
                        out[:] += aux

                    output('eval_einsum NM: {} s'.format(tt.stop()))

            else: # This never happens?
                raise RuntimeError('Impossible code path!')

        output('einsum setup: {} s'.format(timer.stop()))

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

class EDivTerm(ETermBase, Term):
    name = 'dw_ediv'
    arg_types = ('opt_material', 'virtual')
    arg_shapes = [{'opt_material' : '1, 1', 'virtual' : ('D', None)},
                  {'opt_material' : None}]

    def expression(self, mat, virtual, mode=None, term_mode=None,
                   diff_var=None, **kwargs):
        if mat is None:
            expr = self.einsum('i.i', virtual,
                               diff_var=diff_var)

        else:
            expr = self.einsum('0,i.i', mat, virtual,
                               diff_var=diff_var)

        return expr

register_term(EDivTerm)

class EStokesTerm(ETermBase, Term):
    name = 'dw_estokes'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'state', 'virtual'),
                 ('opt_material', 'parameter_v', 'parameter_s'))
    arg_shapes = [{'opt_material' : '1, 1',
                   'virtual/grad' : ('D', None), 'state/grad' : 1,
                   'virtual/div' : (1, None), 'state/div' : 'D',
                   'parameter_v' : 'D', 'parameter_s' : 1},
                  {'opt_material' : None}]
    modes = ('grad', 'div', 'eval')

    def expression(self, coef, vvar, svar, mode=None, term_mode=None,
                   diff_var=None, **kwargs):
        if coef is None:
            if self.mode == 'grad':
                expr = self.einsum('i.i,0', vvar, svar,
                                   diff_var=diff_var)

            else:
                expr = self.einsum('0,i.i', svar, vvar,
                                   diff_var=diff_var)

        else:
            if self.mode == 'grad':
                expr = self.einsum('0,i.i,0', coef, virtual, state,
                                   diff_var=diff_var)

            else:
                expr = self.einsum('0,0,i.i', coef, virtual, state,
                                   diff_var=diff_var)

        return expr

register_term(EStokesTerm)

import numpy as nm

try:
    import dask.array as da

except ImportError:
    da = None

try:
    import opt_einsum as oe

except ImportError:
    oe = None

from pyparsing import (Word, Suppress, oneOf, OneOrMore, delimitedList,
                       Combine, alphas, Literal)

from sfepy.base.base import output, Struct
from sfepy.base.timing import Timer
from sfepy.discrete import FieldVariable
from sfepy.mechanics.tensors import dim2sym
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
    _aux_letters = 'rstuvwxyz'

    def __init__(self, n_add, dc_type, region, dofs_cache):
        self.n_add = n_add
        self.subscripts = [[] for ia in range(n_add)]
        self.operands = [[] for ia in range(n_add)]
        self.operand_names = [[] for ia in range(n_add)]
        self.out_subscripts = ['c' for ia in range(n_add)]
        self.ia = 0
        self.dc_type = dc_type
        self.region = region
        self.dofs_cache = dofs_cache
        self.aux_letters = iter(self._aux_letters)

    @staticmethod
    def make_psg(dim):
        sym = dim2sym(dim)
        psg = nm.zeros((dim, dim, sym))
        if dim == 3:
            psg[0, [0,1,2], [0,3,4]] = 1
            psg[1, [0,1,2], [3,1,5]] = 1
            psg[2, [0,1,2], [4,5,2]] = 1

        elif dim == 2:
            psg[0, [0,1], [0,2]] = 1
            psg[1, [0,1], [2,1]] = 1

        return psg

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

    def add_psg(self, iic, ein, psg, iia=None):
        append_all(self.subscripts, '{}{}{}'.format(iic, ein[2], ein[0]),
                   ii=iia)
        append_all(self.operands, psg, ii=iia)
        append_all(self.operand_names, 'Psg', ii=iia)

    def add_arg_dofs(self, iin, ein, arg, iia=None):
        dofs = self.dofs_cache.get(arg.name)
        if dofs is None:
            conn = arg.field.get_econn(self.dc_type, self.region)
            dofs_vec = arg().reshape((-1, arg.n_components))
            # axis 0: cells, axis 1: node, axis 2: component
            dofs = dofs_vec[conn]
            self.dofs_cache[arg.name] = dofs

        if arg.n_components > 1:
            term = 'c{}{}'.format(iin, ein[0])

        else:
            dofs.shape = (dofs.shape[0], -1)
            term = 'c{}'.format(iin)

        append_all(self.subscripts, term, ii=iia)
        append_all(self.operands, dofs, ii=iia)
        append_all(self.operand_names, arg.name + '.dofs', ii=iia)

    def add_virtual_arg(self, arg, ii, ein, qsb, qsbg, modifier):
        iin = self.letters[ii] # node (qs basis index)
        if ('.' in ein) or (':' in ein): # derivative, symmetric gradient
            self.add_bfg(iin, ein, qsbg, arg.name + '.bfg')

        else:
            self.add_bf(iin, ein, qsb, arg.name + '.bf')

        out_letters = iin

        if arg.n_components > 1:
            iic = next(self.aux_letters) # component
            if ':' not in ein:
                ee = nm.eye(arg.n_components)
                self.add_eye(iic, ein, ee)

            else: # symmetric gradient
                if modifier[0][0] == 's': # vector storage
                    psg = self.make_psg(arg.dim)
                    self.add_psg(iic, ein, psg)

                else:
                    raise ValueError('unknown argument modifier! ({})'
                                     .format(modifier))

            out_letters = iic + out_letters

        for iia in range(self.n_add):
            self.out_subscripts[iia] += out_letters

    def add_state_arg(self, arg, ii, ein, qsb, qsbg, modifier, diff_var):
        iin = self.letters[ii] # node (qs basis index)
        if ('.' in ein) or (':' in ein): # derivative, symmetric gradient
            self.add_bfg(iin, ein, qsbg, arg.name + '.bfg')

        else:
            self.add_bf(iin, ein, qsb, arg.name + '.bf')

        out_letters = iin

        if (diff_var != arg.name):
            if ':' not in ein:
                self.add_arg_dofs(iin, ein, arg)

            else: # symmetric gradient
                if modifier[0][0] == 's': # vector storage
                    iic = next(self.aux_letters) # component
                    psg = self.make_psg(arg.dim)
                    self.add_psg(iic, ein, psg)
                    self.add_arg_dofs(iin, [iic], arg)

                else:
                    raise ValueError('unknown argument modifier! ({})'
                                     .format(modifier))

        else:
            if arg.n_components > 1:
                iic = next(self.aux_letters) # component
                if ':' not in ein:
                    ee = nm.eye(arg.n_components)

                else: # symmetric gradient
                    if modifier[0][0] == 's': # vector storage
                        psg = self.make_psg(arg.dim)

                    else:
                        raise ValueError('unknown argument modifier! ({})'
                                         .format(modifier))

                out_letters = iic + out_letters

            for iia in range(self.n_add):
                if iia != self.ia:
                    self.add_arg_dofs(iin, ein, arg, iia)

                elif arg.n_components > 1:
                    if ':' not in ein:
                        self.add_eye(iic, ein, ee, iia)

                    else:
                        self.add_psg(iic, ein, psg, iia)

            self.out_subscripts[self.ia] += out_letters
            self.ia += 1

    def add_material_arg(self, arg, ii, ein, name):
        append_all(self.subscripts, 'cq{}'.format(ein))
        append_all(self.operands, arg)
        append_all(self.operand_names, name)

    @staticmethod
    def join_subscripts(subscripts, out_subscripts):
        return ','.join(subscripts) + '->' + out_subscripts

    def get_expressions(self):
        expressions = [self.join_subscripts(self.subscripts[ia],
                                            self.out_subscripts[ia])
                       for ia in range(self.n_add)]
        return expressions

    def get_sizes(self, ia):
        return get_sizes(self.subscripts[ia], self.operands[ia])

    def get_output_shape(self, ia):
        return tuple(self.get_sizes(ia)[ii] for ii in self.out_subscripts[ia])

    def print_shapes(self):
        for ia in range(self.n_add):
            sizes = self.get_sizes(ia)
            output(sizes)
            out_shape = self.get_output_shape(ia)
            output(self.out_subscripts[ia], out_shape, '=')

            for name, ii, op in zip(self.operand_names[ia],
                                    self.subscripts[ia],
                                    self.operands[ia]):
                output('  {:10}{:8}{}'.format(name, ii, op.shape))

def collect_modifiers(modifiers):
    def _collect_modifiers(toks):
        if len(toks) > 1:
            out = []
            modifiers.append([])
            for ii, mod in enumerate(toks[::3]):
                tok = toks[3*ii+1]
                tok = tok.replace(tok[0], toks[2])
                modifiers[-1].append(list(toks))
                out.append(tok)
            return out

        else:
            modifiers.append(None)
            return toks
    return _collect_modifiers

def parse_term_expression(texpr):
    mods = 's'
    lparen, rparen = map(Suppress, '()')
    simple_arg = Word(alphas + '.:0')
    arrow = Literal('->').suppress()
    letter = Word(alphas, exact=1)
    mod_arg = oneOf(mods) + lparen + simple_arg + rparen + arrow + letter
    arg = OneOrMore(simple_arg ^ mod_arg)
    modifiers = []
    arg.setParseAction(collect_modifiers(modifiers))

    parser = delimitedList(Combine(arg))
    eins = parser.parseString(texpr, parseAll=True)
    return eins, modifiers

class ETermBase(Struct):
    """
    Reserved letters:

    c .. cells
    q .. quadrature points
    d-h .. DOFs axes
    r-z .. auxiliary axes
    """
    verbose = 0

    can_backend = {
        'numpy' : nm,
        'numpy_loop' : nm,
        'opt_einsum' : oe,
        'opt_einsum_loop' : oe,
        'dask_single' : da,
        'dask_threads' : da,
        'opt_einsum_dask_single' : oe and da,
        'opt_einsum_dask_threads' : oe and da,
    }

    def set_backend(self, backend='numpy', optimize=True):
        if backend not in self.can_backend.keys():
            raise ValueError('backend {} not in {}!'
                             .format(self.backend, self.can_backend.keys()))

        if not self.can_backend[backend]:
            raise ValueError('backend {} is not available!'.format(backend))

        self.backend = backend
        self.optimize = optimize
        self.paths, self.path_infos = None, None

    def build_expression(self, texpr, *args, diff_var=None):
        timer = Timer('')
        timer.start()

        if diff_var is not None:
            n_add = len([arg.name for arg in args
                         if (isinstance(arg, FieldVariable)
                             and (arg.name == diff_var))])

        else:
            n_add = 1

        eins, modifiers = parse_term_expression(texpr)

        dofs_cache = {}
        self.ebuilder = ExpressionBuilder(
            n_add, self.get_dof_conn_type(), self.region, dofs_cache,
        )

        # Virtual variable must be the first variable.
        # Numpy arrays cannot be compared -> use a loop.
        for iv, arg in enumerate(args):
            if isinstance(arg, FieldVariable) and arg.is_virtual():
                ag, _ = self.get_mapping(arg)
                self.ebuilder.add_constant(ag.det[..., 0, 0], 'J')
                self.ebuilder.add_virtual_arg(arg, iv, eins[iv], ag.bf,
                                              ag.bfg, modifiers[iv])
                break
        else:
            iv = -1
            for ip, arg in enumerate(args):
                if (isinstance(arg, FieldVariable)
                    and arg.is_state_or_parameter()):
                    ag, _ = self.get_mapping(arg)
                    self.ebuilder.add_constant(ag.det[..., 0, 0], 'J')
                    break
            else:
                raise ValueError('no FieldVariable in arguments!')

        for ii, ein in enumerate(eins):
            if ii == iv: continue
            arg = args[ii]

            if isinstance(arg, nm.ndarray):
                self.ebuilder.add_material_arg(arg, ii, ein,
                                               '.'.join(self.arg_names[ii]))

            elif isinstance(arg, FieldVariable) and arg.is_state():
                ag, _ = self.get_mapping(arg)
                self.ebuilder.add_state_arg(arg, ii, ein, ag.bf, ag.bfg,
                                            modifiers[ii], diff_var)

            else:
                raise ValueError('unknown argument type! ({})'
                                 .format(type(arg)))

        if self.verbose:
            output('build expression: {} s'.format(timer.stop()))

    def get_paths(self, expressions, operands):
        if ('numpy' in self.backend) or self.backend.startswith('dask'):
            paths, path_infos = zip(*[nm.einsum_path(
                expressions[ia], *operands[ia],
                optimize=self.optimize,
            ) for ia in range(len(operands))])

        elif 'opt_einsum' in self.backend:
            paths, path_infos = zip(*[oe.contract_path(
                expressions[ia], *operands[ia],
                optimize=self.optimize,
            ) for ia in range(len(operands))])

        else:
            raise ValueError('unsupported backend! ({})'.format(self.backend))

        return paths, path_infos

    def make_function(self, texpr, *args, diff_var=None):
        timer = Timer('')
        timer.start()

        if not hasattr(self, 'ebuilder'):
            self.build_expression(texpr, *args, diff_var=diff_var)

        if not hasattr(self, 'paths') or (self.paths is None):
            self.parsed_expressions = self.ebuilder.get_expressions()
            if self.verbose:
                output(self.parsed_expressions)
            if self.verbose > 1:
                self.ebuilder.print_shapes()

            self.paths, self.path_infos = self.get_paths(
                self.parsed_expressions,
                self.ebuilder.operands,
            )
            if self.verbose > 2:
                for path, path_info in zip(self.paths, self.path_infos):
                    output(path)
                    output(path_info)

        operands = self.ebuilder.operands
        n_add = len(operands)

        if 'numpy' in self.backend:
            def eval_einsum(out, eshape):
                vout = out.reshape(eshape)
                nm.einsum(self.parsed_expressions[0], *operands[0],
                          out=vout,
                          optimize=self.paths[0])
                for ia in range(1, n_add):
                    aux = nm.einsum(self.parsed_expressions[ia],
                                    *operands[ia],
                                    optimize=self.paths[ia])
                    out[:] += aux.reshape(out.shape)

        elif 'opt_einsum' in self.backend:
            def eval_einsum(out, eshape):
                vout = out.reshape(eshape)
                oe.contract(self.parsed_expressions[0], *operands[0],
                            out=vout,
                            optimize=self.paths[0])
                for ia in range(1, n_add):
                    aux = oe.contract(self.parsed_expressions[ia],
                                      *operands[ia],
                                      optimize=self.paths[ia])
                    out[:] += aux.reshape(out.shape)

        elif self.backend.startswith('dask'):
            scheduler = {'dask_single' : 'single-threaded',
                         'dask_threads' : 'threads'}[self.backend]
            def eval_einsum(out, eshape):
                vout = out.reshape(eshape)
                _out = da.einsum(self.parsed_expressions[0], *operands[0],
                                 out=vout,
                                 optimize=self.paths[0])
                for ia in range(1, n_add):
                    aux = da.einsum(self.parsed_expressions[ia],
                                    *operands[ia],
                                    optimize=self.paths[ia])
                    _out += aux

                _out.compute(scheduler=scheduler)

        else:
            raise ValueError('unsupported backend! ({})'.format(self.backend))

        if self.verbose:
            output('einsum setup: {} s'.format(timer.stop()))

        return eval_einsum

    @staticmethod
    def function(out, eval_einsum, eshape):
        tt = Timer('')
        tt.start()
        eval_einsum(out, eshape)
        output('eval_einsum: {} s'.format(tt.stop()))
        return 0

    def get_fargs(self, *args, **kwargs):
        mode, term_mode, diff_var = args[-3:]

        if mode == 'weak':
            vvar = self.get_virtual_variable()
            n_elr, n_qpr, dim, n_enr, n_cr = self.get_data_shape(vvar)

            if diff_var is not None:
                varc = self.get_variables(as_list=False)[diff_var]
                n_elc, n_qpc, dim, n_enc, n_cc = self.get_data_shape(varc)
                eshape = tuple([n_elr]
                               + ([n_cr] if n_cr > 1 else [])
                               + [n_enr]
                               + ([n_cc] if n_cc > 1 else [])
                               + [n_enc])

            else:
                eshape = (n_elr, n_cr, n_enr) if n_cr > 1 else (n_elr, n_enr)

        else:
            if diff_var is not None:
                raise ValueError('cannot differentiate in {} mode!'
                                 .format(mode))

            # self.ebuilder is created in self.get_eval_shape() call by Term.
            eshape = self.ebuilder.get_output_shape(0)

        eval_einsum = self.get_function(*args, **kwargs)

        return eval_einsum, eshape

    def get_eval_shape(self, *args, **kwargs):
        mode, term_mode, diff_var = args[-3:]
        if diff_var is not None:
            raise ValueError('cannot differentiate in {} mode!'
                             .format(mode))

        self.get_function(*args, **kwargs)

        out_shape = self.ebuilder.get_output_shape(0)

        operands = self.ebuilder.operands[0]
        dtype = nm.find_common_type([op.dtype for op in operands], [])

        return out_shape, dtype

class ELaplaceTerm(ETermBase, Term):
    name = 'dw_elaplace'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'parameter_1', 'parameter_2'))
    arg_shapes = [{'opt_material' : '1, 1', 'virtual' : (1, 'state'),
                   'state' : 1, 'parameter_1' : 1, 'parameter_2' : 1},
                  {'opt_material' : None}]
    modes = ('weak', 'eval')

    def get_function(self, mat, virtual, state, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        """
        diff_var not needed here(?), but Term passes it in *args.
        """
        if mat is None:
            fun = self.make_function(
                '0.j,0.j', virtual, state, diff_var=diff_var,
            )

        else:
            fun = self.make_function(
                'jk,0.j,0.k', mat, virtual, state, diff_var=diff_var,
            )

        return fun

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

    def get_function(self, mat, virtual, state, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        if mat is None:
            fun = self.make_function(
                'i,i', virtual, state, diff_var=diff_var,
            )

        else:
            fun = self.make_function(
                'ij,i,j', mat, virtual, state, diff_var=diff_var,
            )

        return fun

register_term(EVolumeDotTerm)

class EConvectTerm(ETermBase, Term):
    name = 'dw_econvect'
    arg_types = ('virtual', 'state')
    arg_shapes = {'virtual' : ('D', 'state'), 'state' : 'D'}

    def get_function(self, virtual, state, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        return self.make_function(
            'i,i.j,j', virtual, state, state, diff_var=diff_var,
        )

register_term(EConvectTerm)

class EDivTerm(ETermBase, Term):
    name = 'dw_ediv'
    arg_types = ('opt_material', 'virtual')
    arg_shapes = [{'opt_material' : '1, 1', 'virtual' : ('D', None)},
                  {'opt_material' : None}]

    def get_function(self, mat, virtual, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        if mat is None:
            fun = self.make_function(
                'i.i', virtual, diff_var=diff_var,
            )

        else:
            fun = self.make_function(
                '0,i.i', mat, virtual, diff_var=diff_var,
            )

        return fun

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

    def get_function(self, coef, var_v, var_s, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        if coef is None:
            fun = self.make_function(
                'i.i,0', var_v, var_s, diff_var=diff_var,
            )

        else:
            fun = self.make_function(
                '0,i.i,0', coef, var_v, var_s, diff_var=diff_var,
            )

        return fun

register_term(EStokesTerm)

class ELinearElasticTerm(ETermBase, Term):
    name = 'dw_elin_elastic'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'parameter_1', 'parameter_2'))
    arg_shapes = {'material' : 'S, S', 'virtual' : ('D', 'state'),
                  'state' : 'D', 'parameter_1' : 'D', 'parameter_2' : 'D'}
    modes = ('weak', 'eval')

    def get_function(self, mat, virtual, state, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        return self.make_function(
            'IK,s(i:j)->I,s(k:l)->K', mat, virtual, state, diff_var=diff_var,
        )

register_term(ELinearElasticTerm)

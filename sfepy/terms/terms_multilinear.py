import numpy as nm
from sfepy.linalg import dot_sequences

try:
    import dask.array as da

except ImportError:
    da = None

try:
    import opt_einsum as oe

except ImportError:
    oe = None

try:
    from jax import config
    config.update("jax_enable_x64", True)
    import jax
    import jax.numpy as jnp

except ImportError:
    jnp = jax = None

from pyparsing import (Word, Suppress, oneOf, OneOrMore, delimitedList,
                       Combine, alphas, alphanums, Literal)
from functools import partial
from copy import copy

from sfepy.base.base import output, Struct
from sfepy.base.timing import Timer
from sfepy.mechanics.tensors import dim2sym
from sfepy.terms.terms import Term

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
    mods = 's', 'v'
    lparen, rparen = map(Suppress, '()')
    simple_arg = Word(alphanums + '.:')
    arrow = Literal('->').suppress()
    letter = Word(alphas, exact=1)
    mod_arg = oneOf(mods) + lparen + simple_arg + rparen + arrow + letter
    arg = OneOrMore(simple_arg ^ mod_arg)
    modifiers = []
    arg.setParseAction(collect_modifiers(modifiers))

    parser = delimitedList(Combine(arg))
    eins = parser.parseString(texpr, parseAll=True)
    return eins, modifiers

def append_all(seqs, item, ii=None):
    if ii is None:
        for seq in seqs:
            seq.append(copy(item))

    else:
        seqs[ii].append(copy(item))

def get_sizes(indices, operands):
    sizes = {}

    for iis, op in zip(indices, operands):
        for ii, size in zip(iis, op.shape):
            sizes[ii] = size

    return sizes

def get_output_shape(out_subscripts, subscripts, operands):
    return tuple(get_sizes(subscripts, operands)[ii] for ii in out_subscripts)

def find_free_indices(indices):
    ii = ''.join(indices)
    ifree = [c for c in set(ii) if ii.count(c) == 1]
    return ifree

def get_loop_indices(subs, loop_index):
    return [indices.index(loop_index) if loop_index in indices else None
            for indices in subs]

def get_einsum_ops(eargs, ebuilder, expr_cache):
    dargs = {arg.name : arg for arg in eargs}

    operands = [[] for ia in range(ebuilder.n_add)]
    for ia in range(ebuilder.n_add):
        for io, oname in enumerate(ebuilder.operand_names[ia]):
            arg_name, val_name = oname.split('.')
            arg = dargs[arg_name]
            if val_name == 'dofs':
                step_cache = arg.arg.evaluate_cache.setdefault('dofs', {})
                step = arg.term.arg_steps[arg.name]
                cache = step_cache.setdefault(step, {})
                op = arg.get_dofs(cache, expr_cache, oname)

            elif val_name == 'bf':
                cell_dep = len(ebuilder.subscripts[ia][io]) == 3
                op = arg.get_bf(expr_cache, cell_dependent=cell_dep)

            elif val_name == 'bfg':
                ag, _ = arg.term.get_mapping(arg.arg)
                op = ag.bfg

            elif val_name == 'det':
                ag, _ = arg.term.get_mapping(arg.arg)
                op = ag.det[..., 0, 0]

            elif val_name == 'ivol':
                ag, _ = arg.term.get_mapping(arg.arg)
                op = 1.0 / ag.volume[:, 0, 0, 0]

            elif val_name == 'I':
                op = ebuilder.make_eye(arg.n_components)

            elif val_name == 'Psg':
                op = ebuilder.make_psg(arg.dim)

            elif val_name == 'Pvg':
                op = ebuilder.make_pvg(arg.dim)

            else:
                op = dargs[arg_name].get(
                    val_name,
                    msg_if_none='{} has no attribute {}!'
                    .format(arg_name, val_name)
                )
                ics = ebuilder.components[ia][io]
                if len(ics):
                    iis = [slice(None)] * 2
                    iis += [slice(None) if ic is None else ic for ic in ics]
                    op = op[tuple(iis)]

            operands[ia].append(op)

    return operands

def get_slice_ops(subs, ops, loop_index):
    ics = get_loop_indices(subs, loop_index)

    def slice_ops(ic):
        sops = []
        for ii, icol in enumerate(ics):
            op = ops[ii]
            if icol is not None:
                slices = tuple(slice(None, None) if isub != icol else ic
                               for isub in range(op.ndim))
                sops.append(op[slices])

            else:
                sops.append(op)

        return sops

    return slice_ops

class ExpressionArg(Struct):

    @staticmethod
    def from_term_arg(arg, term):
        from sfepy.discrete import FieldVariable

        if isinstance(arg, ExpressionArg):
            return arg

        if isinstance(arg, FieldVariable):
            obj = ExpressionArg(name=arg.name, arg=arg, term=term,
                                n_components=arg.n_components, dim=arg.dim,
                                kind='virtual' if arg.is_virtual() else 'state')

        elif isinstance(arg, nm.ndarray):
            aux = term.get_args()
            # Find arg in term arguments using a loop (numpy arrays cannot be
            # compared) to get its name.
            ii = [ii for ii in range(len(term.args)) if aux[ii] is arg][0]
            obj = ExpressionArg(name='_'.join(term.arg_names[ii]), arg=arg,
                                kind='ndarray')

        elif isinstance(arg, tuple) and isinstance(arg[0], nm.ndarray):
            obj = ExpressionArg(name=arg[1], arg=arg[0], kind='ndarray')

        else:
            raise ValueError('unknown argument type! ({})'.format(type(arg)))

        return obj

    def get_dofs(self, cache, expr_cache, oname):
        if self.kind != 'state': return

        key = (self.name, self.term.region.name,
               self.term.arg_derivatives[self.name])
        dofs = cache.get(key)
        if dofs is None:
            arg = self.arg
            dofs_vec = self.arg(derivative=key[-1])
            dofs_vec = dofs_vec.reshape((-1, arg.n_components))
            # # axis 0: cells, axis 1: node, axis 2: component
            # dofs = dofs_vec[conn]
            # axis 0: cells, axis 1: component, axis 2: node
            conn = arg.field.get_econn(self.term.get_dof_conn_type(arg.name),
                                       self.term.region)
            dofs = dofs_vec[conn].transpose((0, 2, 1))
            if arg.n_components == 1:
                dofs.shape = (dofs.shape[0], -1)
            cache[key] = dofs

            # New dofs -> clear dofs from expression cache.
            for key in list(expr_cache.keys()):
                if isinstance(key, tuple) and (key[0] == oname):
                    expr_cache.pop(key)

        return dofs

    def get_bf(self, expr_cache, cell_dependent=False):
        ag, _ = self.term.get_mapping(self.arg)
        if self.term.integration == 'facet_extra':
            if 'L2_constant' in self.arg.field.family_name:
                bf = ag.bf[:, :1]

            else:
                sd = self.arg.field.extra_data[f'sd_{self.term.region.name}']
                _bf = self.arg.field.eval_basis(sd.bkey, 0, self.term.integral)
                bf = _bf[sd.fis[:, 1], ...]

        else:
            bf = ag.bf

        key = 'bf{}'.format(id(bf))
        _bf  = expr_cache.get(key)

        if cell_dependent:
            if _bf is None:
                _bf = bf[:, :, 0]
                expr_cache[key] = _bf

        else:
            if _bf is None:
                _bf = bf[0, :, 0]
                expr_cache[key] = _bf

        return _bf

class ExpressionBuilder(Struct):
    letters = 'defgh'
    _aux_letters = 'rstuvwxyz'

    def __init__(self, n_add, cache):
        self.n_add = n_add
        self.subscripts = [[] for ia in range(n_add)]
        self.operand_names = [[] for ia in range(n_add)]
        self.components = [[] for ia in range(n_add)]
        self.out_subscripts = ['c' for ia in range(n_add)]
        self.ia = 0
        self.cache = cache
        self.aux_letters = iter(self._aux_letters)

    def make_eye(self, size):
        key = 'I{}'.format(size)
        ee = self.cache.get(key)
        if ee is None:
            ee = nm.eye(size)
            self.cache[key] = ee

        return ee

    def make_psg(self, dim):
        key = 'Psg{}'.format(dim)
        psg = self.cache.get(key)
        if psg is None:
            sym = dim2sym(dim)
            psg = nm.zeros((dim, dim, sym))
            if dim == 3:
                psg[0, [0,1,2], [0,3,4]] = 1
                psg[1, [0,1,2], [3,1,5]] = 1
                psg[2, [0,1,2], [4,5,2]] = 1

            elif dim == 2:
                psg[0, [0,1], [0,2]] = 1
                psg[1, [0,1], [2,1]] = 1

            self.cache[key] = psg

        return psg

    def make_pvg(self, dim):
        key = 'Pvg{}'.format(dim)
        pvg = self.cache.get(key)
        if pvg is None:
            pvg = nm.zeros((dim, dim, dim*dim))
            if dim == 3:
                pvg[0, [0, 1, 2], [0, 1, 2]] = 1
                pvg[1, [0, 1, 2], [3, 4, 5]] = 1
                pvg[2, [0, 1, 2], [6, 7, 8]] = 1

            elif dim == 2:
                pvg[0, [0, 1], [0, 1]] = 1
                pvg[1, [0, 1], [2, 3]] = 1

            self.cache[key] = pvg

        return pvg

    def add_cell_scalar(self, name, cname):
        append_all(self.subscripts, 'c')
        append_all(self.operand_names, '.'.join((name, cname)))
        append_all(self.components, [])

    def add_qp_scalar(self, name, cname):
        append_all(self.subscripts, 'cq')
        append_all(self.operand_names, '.'.join((name, cname)))
        append_all(self.components, [])

    def add_bfg(self, iin, ein, name):
        append_all(self.subscripts, 'cq{}{}'.format(ein[2], iin))
        append_all(self.operand_names, name + '.bfg')
        append_all(self.components, [])

    def add_bf(self, iin, ein, name, cell_dependent=False):
        if cell_dependent:
            append_all(self.subscripts, 'cq{}'.format(iin))

        else:
            append_all(self.subscripts, 'q{}'.format(iin))

        append_all(self.operand_names, name + '.bf')
        append_all(self.components, [])

    def add_eye(self, iic, ein, name, iia=None):
        append_all(self.subscripts, '{}{}'.format(ein[0], iic), ii=iia)
        append_all(self.operand_names, '{}.I'.format(name), ii=iia)
        append_all(self.components, [])

    def add_psg(self, iic, ein, name, iia=None):
        append_all(self.subscripts, '{}{}{}'.format(iic, ein[2], ein[0]),
                   ii=iia)
        append_all(self.operand_names, name + '.Psg', ii=iia)
        append_all(self.components, [])

    def add_pvg(self, iic, ein, name, iia=None):
        append_all(self.subscripts, '{}{}{}'.format(iic, ein[2], ein[0]),
                   ii=iia)
        append_all(self.operand_names, name + '.Pvg', ii=iia)
        append_all(self.components, [])

    def add_arg_dofs(self, iin, ein, name, n_components, iia=None):
        if n_components > 1:
            #term = 'c{}{}'.format(iin, ein[0])
            term = 'c{}{}'.format(ein[0], iin)

        else:
            term = 'c{}'.format(iin)

        append_all(self.subscripts, term, ii=iia)
        append_all(self.operand_names, name + '.dofs', ii=iia)
        append_all(self.components, [])

    def add_virtual_arg(self, arg, ii, ein, modifier):
        iin = self.letters[ii] # node (qs basis index)
        if ('.' in ein) or (':' in ein): # derivative, symmetric gradient
            self.add_bfg(iin, ein, arg.name)

        else:
            cell_dep = arg.term.integration == 'facet_extra'
            self.add_bf(iin, ein, arg.name, cell_dependent=cell_dep)

        out_letters = iin

        if arg.n_components > 1:
            iic = next(self.aux_letters) # component

            if modifier is not None:
                # symmetric gradient, vector storage
                if ':' in ein and modifier[0][0] == 's':
                    self.add_psg(iic, ein, arg.name)
                # nonsymmetric gradient, vector storage
                elif modifier[0][0] == 'v':
                    self.add_pvg(iic, ein, arg.name)
                else:
                    raise ValueError('unknown argument modifier! ({})'
                                     .format(modifier))
            else:
                self.add_eye(iic, ein, arg.name)

            out_letters = iic + out_letters

        for iia in range(self.n_add):
            self.out_subscripts[iia] += out_letters

    def add_state_arg(self, arg, ii, ein, modifier, diff_var):
        iin = self.letters[ii] # node (qs basis index)
        if ('.' in ein) or (':' in ein): # derivative, symmetric gradient
            self.add_bfg(iin, ein, arg.name)

        else:
            cell_dep = arg.term.integration == 'facet_extra'
            self.add_bf(iin, ein, arg.name, cell_dependent=cell_dep)

        out_letters = iin

        if (diff_var != arg.name):
            if modifier is not None:
                # symmetric gradient, vector storage
                if ':' in ein and modifier[0][0] == 's':
                    iic = next(self.aux_letters)  # component
                    self.add_psg(iic, ein, arg.name)
                    self.add_arg_dofs(iin, [iic], arg.name, arg.n_components)
                # nonsymmetric gradient, vector storage
                elif modifier[0][0] == 'v':
                    iic = next(self.aux_letters)  # component
                    self.add_pvg(iic, ein, arg.name)
                    self.add_arg_dofs(iin, [iic], arg.name, arg.n_components)
                else:
                    raise ValueError('unknown argument modifier! ({})'
                                     .format(modifier))
            else:
                self.add_arg_dofs(iin, ein, arg.name, arg.n_components)

        else:
            if arg.n_components > 1:
                iic = next(self.aux_letters)  # component
                out_letters = iic + out_letters

            for iia in range(self.n_add):
                if iia != self.ia:
                    self.add_arg_dofs(iin, ein, arg.name, arg.n_components, iia)

                elif arg.n_components > 1:
                    if modifier is not None:
                        if ':' in ein and modifier[0][0] == 's':
                            self.add_psg(iic, ein, arg.name, iia)
                        elif modifier[0][0] == 'v':
                            self.add_pvg(iic, ein, arg.name, iia)
                        else:
                            raise ValueError('unknown argument modifier! ({})'
                                             .format(modifier))
                    else:
                        self.add_eye(iic, ein, arg.name, iia)

            self.out_subscripts[self.ia] += out_letters
            self.ia += 1

    def add_material_arg(self, arg, ii, ein):
        append_all(self.components, [])
        rein = []
        for ii, ie in enumerate(ein):
            if str.isnumeric(ie):
                for comp in self.components:
                    comp[-1].append(int(ie))

            else:
                for comp in self.components:
                    comp[-1].append(None)
                rein.append(ie)
        rein = ''.join(rein)

        append_all(self.subscripts, 'cq{}'.format(rein))
        append_all(self.operand_names, arg.name + '.arg')

    def build(self, texpr, *args, mode=None, diff_var=None):
        eins, modifiers = parse_term_expression(texpr)

        # Virtual variable must be the first variable.
        # Numpy arrays cannot be compared -> use a loop.
        for iv, arg in enumerate(args):
            if arg.kind == 'virtual':
                if mode != 'qp':
                    self.add_qp_scalar(arg.name, 'det')
                self.add_virtual_arg(arg, iv, eins[iv], modifiers[iv])
                break
        else:
            iv = -1
            for ip, arg in enumerate(args):
                if arg.kind == 'state':
                    if mode != 'qp':
                        self.add_qp_scalar(arg.name, 'det')
                        if mode == 'el_avg':
                            self.add_cell_scalar(arg.name, 'ivol')
                    break
            else:
                raise ValueError('no FieldVariable in arguments!')

        for ii, ein in enumerate(eins):
            if ii == iv: continue
            arg = args[ii]

            if arg.kind == 'ndarray':
                self.add_material_arg(arg, ii, ein)

            elif arg.kind == 'state':
                self.add_state_arg(arg, ii, ein, modifiers[ii], diff_var)

            else:
                raise ValueError('unknown argument type! ({})'
                                 .format(type(arg)))

        for ia, subscripts in enumerate(self.subscripts):
            ifree = [ii for ii in find_free_indices(subscripts)
                     if ii not in self.out_subscripts[ia]]
            if ifree:
                # Lexicographic ordering of output indices, i.e.
                # (n_comp, dim) or (dim, n_comp) <=> i.j or j.i
                self.out_subscripts[ia] += ''.join(sorted(ifree))

            if mode == 'qp':
                self.out_subscripts[ia] = (self.out_subscripts[ia]
                                           .replace('c', 'cq'))

    @staticmethod
    def join_subscripts(subscripts, out_subscripts):
        return ','.join(subscripts) + '->' + out_subscripts

    def get_expressions(self, subscripts=None):
        if subscripts is None:
            subscripts = self.subscripts

        expressions = [self.join_subscripts(subscripts[ia],
                                            self.out_subscripts[ia])
                       for ia in range(self.n_add)]
        return tuple(expressions)

    def print_shapes(self, subscripts, operands):
        if subscripts is None:
            subscripts = self.subscripts
        output('number of expressions:', self.n_add)
        for onames, outs, subs, ops in zip(
            self.operand_names, self.out_subscripts, subscripts, operands,
        ):
            sizes = get_sizes(subs, ops)
            output(sizes)
            out_shape = get_output_shape(outs, subs, ops)
            output(outs, out_shape, '=')

            for name, ii, op in zip(onames, subs, ops):
                output('  {:10} {:8} {}'.format(name, ii, op.shape))

    def apply_layout(self, layout, operands, defaults=None, verbosity=0):
        if layout == 'cqgvd0':
            return self.subscripts, operands

        if defaults is None:
            defaults = {
                'det' : 'cq',
                'bf' : ('qd', 'cqd'),
                'bfg' : 'cqgd',
                'dofs' : ('cd', 'cvd'),
                'mat' : 'cq',
            }

        mat_range = ''.join([str(ii) for ii in range(10)])
        new_subscripts = [subs.copy() for subs in self.subscripts]
        new_operands = [ops.copy() for ops in operands]
        for ia in range(self.n_add):
            for io, (oname, subs, op) in enumerate(zip(self.operand_names[ia],
                                                      self.subscripts[ia],
                                                      operands[ia])):
                arg_name, val_name = oname.split('.')

                if val_name in ('det','bfg'):
                    default = defaults[val_name]

                elif val_name in ('bf', 'dofs'):
                    default = defaults[val_name][op.ndim - 2]

                elif val_name in ('I', 'Psg', 'Pvg'):
                    default = layout.replace('0', '') # -> Do nothing.

                else:
                    default = defaults['mat'] + mat_range[:(len(subs) - 2)]

                if '0' in default: # Material
                    inew = nm.array([default.find(il)
                                     for il in layout.replace('0', default[2:])
                                     if il in default])

                else:
                    inew = nm.array([default.find(il)
                                     for il in layout if il in default])

                new = ''.join([default[ii] for ii in inew])
                if verbosity > 2:
                    output(arg_name, val_name, subs, default, op.shape, layout)
                    output(inew, new)

                if new == default:
                    new_subscripts[ia][io] = subs
                    new_operands[ia][io] = op

                else:
                    new_subs = ''.join([subs[ii] for ii in inew])
                    if val_name == 'dofs':
                        key = (oname,) + tuple(inew)

                    else:
                        # id is unique only during object lifetime!
                        key = (id(op),) + tuple(inew)

                    new_op = self.cache.get(key)
                    if new_op is None:
                        new_op = op.transpose(inew).copy()
                        self.cache[key] = new_op

                    new_subscripts[ia][io] = new_subs
                    new_operands[ia][io] = new_op

                if verbosity > 2:
                    output('->', new_subscripts[ia][io])

        return new_subscripts, new_operands

    def transform(self, subscripts, operands, transformation='loop', **kwargs):
        if transformation == 'loop':
            expressions, poperands, all_slice_ops, loop_sizes = [], [], [], []
            loop_index = kwargs.get('loop_index', 'c')

            for ia, (subs, out_subscripts, ops) in enumerate(zip(
                    subscripts, self.out_subscripts, operands
            )):
                slice_ops = get_slice_ops(subs, ops, loop_index)
                tsubs = [ii.replace(loop_index, '') for ii in subs]
                tout_subs = out_subscripts.replace(loop_index, '')
                expr = self.join_subscripts(tsubs, tout_subs)
                pops = slice_ops(0)
                expressions.append(expr)
                poperands.append(pops)
                all_slice_ops.append(slice_ops)
                loop_sizes.append(get_sizes(subs, ops)[loop_index])

            return tuple(expressions), poperands, all_slice_ops, loop_sizes

        elif transformation == 'dask':
            da_operands = []
            c_chunk_size = kwargs.get('c_chunk_size')
            loop_index = kwargs.get('loop_index', 'c')
            for ia in range(len(operands)):
                da_ops = []
                for name, ii, op in zip(self.operand_names[ia],
                                        subscripts[ia],
                                        operands[ia]):
                    if loop_index in ii:
                        if c_chunk_size is None:
                            chunks = 'auto'

                        else:
                            ic = ii.index(loop_index)
                            chunks = (op.shape[:ic]
                                      + (c_chunk_size,)
                                      + op.shape[ic + 1:])

                        da_op = da.from_array(op, chunks=chunks, name=name)

                    else:
                        da_op = op

                    da_ops.append(da_op)
                da_operands.append(da_ops)

            return da_operands

        else:
            raise ValueError('unknown transformation! ({})'
                             .format(transformation))

class ETermBase(Term):
    """
    Reserved letters:

    c .. cells
    q .. quadrature points
    d-h .. DOFs axes
    r-z .. auxiliary axes

    Layout specification letters:

    c .. cells
    q .. quadrature points
    v .. variable component - matrix form (v, d) -> vector v*d
    g .. gradient component
    d .. local DOF (basis, node)
    0 .. all material axes
    """
    verbosity = 0

    can_backend = {
        'numpy' : nm,
        'numpy_loop' : nm,
        'numpy_qloop' : nm,
        'opt_einsum' : oe,
        'opt_einsum_loop' : oe,
        'opt_einsum_qloop' : oe,
        'jax' : jnp,
        'jax_vmap' : jnp,
        'dask_single' : da,
        'dask_threads' : da,
        'opt_einsum_dask_single' : oe and da,
        'opt_einsum_dask_threads' : oe and da,
    }

    layout_letters = 'cqgvd0'

    def __init__(self, *args, **kwargs):
        Term.__init__(self, *args, **kwargs)

        self.set_verbosity(kwargs.get('verbosity', 0))
        self.set_backend(**kwargs)

    def _eval(self, out, shape, fargs, mode='eval', term_mode=None,
              diff_var=None, **kwargs):
        status = self.function(out, *fargs)
        if mode == 'eval':
            # Sum over elements but not over components.
            out1 = nm.sum(out, 0).squeeze()
            return out1, status

        else:
            return out, status

    def eval_real(self, shape, fargs, mode='eval', term_mode=None,
                  diff_var=None, **kwargs):
        out = nm.empty(shape, dtype=nm.float64)
        return self._eval(out, shape, fargs, mode=mode,
                          term_mode=term_mode, diff_var=diff_var, **kwargs)

    def eval_complex(self, shape, fargs, mode='eval', term_mode=None,
                     diff_var=None, **kwargs):
        out = nm.empty(shape, dtype=nm.complex128)
        return self._eval(out, shape, fargs, mode=mode,
                          term_mode=term_mode, diff_var=diff_var, **kwargs)

    @staticmethod
    def function_timer(out, eval_einsum, *args):
        tt = Timer('')
        tt.start()
        eval_einsum(out, *args)
        output('eval_einsum: {} s'.format(tt.stop()))
        return 0

    @staticmethod
    def function_silent(out, eval_einsum, *args):
        eval_einsum(out, *args)
        return 0

    def set_verbosity(self, verbosity=None):
        if verbosity is not None:
            self.verbosity = verbosity

        if self.verbosity > 0:
            self.function = self.function_timer

        else:
            self.function = self.function_silent

    def set_backend(self, backend='numpy', optimize=True, layout=None,
                    **kwargs):
        if backend not in self.can_backend.keys():
            raise ValueError('backend {} not in {}!'
                             .format(self.backend, self.can_backend.keys()))

        if not self.can_backend[backend]:
            raise ValueError('backend {} is not available!'.format(backend))

        if (hasattr(self, 'backend')
            and (backend == self.backend) and (optimize == self.optimize)
            and (layout == self.layout) and (kwargs == self.backend_kwargs)):
            return

        if layout is not None:
            if set(layout) != set(self.layout_letters):
                raise ValueError('layout can contain only "{}" letters! ({})'
                                 .format(self.layout_letters, layout))
            self.layout = layout

        else:
            self.layout = self.layout_letters

        self.backend = backend
        self.optimize = optimize
        self.backend_kwargs = kwargs
        self.einfos = {}
        self.clear_cache()

    def clear_cache(self):
        self.expr_cache = {}

    def build_expression(self, texpr, *eargs, mode=None, diff_var=None):
        timer = Timer('')
        timer.start()

        if diff_var is not None:
            n_add = len([arg.name for arg in eargs
                         if (arg.kind == 'state')
                         and (arg.name == diff_var)])

        else:
            n_add = 1

        ebuilder = ExpressionBuilder(n_add, self.expr_cache)
        ebuilder.build(texpr, *eargs, mode=mode, diff_var=diff_var)
        if self.verbosity:
            output('build expression: {} s'.format(timer.stop()))

        return ebuilder

    def make_function(self, texpr, *args, mode=None, diff_var=None):
        timer = Timer('')
        timer.start()

        einfo = self.einfos.setdefault(diff_var, Struct(
            eargs=None,
            ebuilder=None,
            paths=None,
            path_infos=None,
            eval_einsum=None,
        ))

        if einfo.eval_einsum is not None:
            if self.verbosity:
                output('einsum setup: {} s'.format(timer.stop()))

            # Use actual material (ndarray) arguments.
            for ii, (arg, earg) in enumerate(zip(args, einfo.eargs)):
                if earg.kind == 'ndarray':
                    if isinstance(arg, tuple):
                        earg.arg = arg[0]
                        earg.name = arg[1]

                    elif isinstance(arg, ExpressionArg):
                        earg.arg = arg.arg

                    else:
                        earg.arg = arg

            return einfo.eval_einsum

        if einfo.eargs is None:
            einfo.eargs = [
                ExpressionArg.from_term_arg(arg, self)
                for arg in args
            ]

        if einfo.ebuilder is None:
            einfo.ebuilder = self.build_expression(texpr, *einfo.eargs,
                                                   mode=mode,
                                                   diff_var=diff_var)

        n_add = einfo.ebuilder.n_add

        if self.backend in ('numpy', 'opt_einsum'):
            contract = nm.einsum if self.backend == 'numpy' else oe.contract
            def eval_einsum_orig(out, eshape, expressions, operands, paths):
                if operands[0][0].flags.c_contiguous:
                    # This is very slow if vout layout differs from operands
                    # layout.
                    vout = out.reshape(eshape)
                    contract(expressions[0], *operands[0], out=vout,
                             optimize=paths[0])

                else:
                    aux = contract(expressions[0], *operands[0],
                                   optimize=paths[0])
                    out[:] += aux.reshape(out.shape)

                for ia in range(1, n_add):
                    aux = contract(expressions[ia], *operands[ia],
                                   optimize=paths[ia])
                    out[:] += aux.reshape(out.shape)

            def eval_einsum0(out, eshape, expressions, operands, paths):
                aux = contract(expressions[0], *operands[0],
                               optimize=paths[0])
                out[:] = aux.reshape(out.shape)
                for ia in range(1, n_add):
                    aux = contract(expressions[ia], *operands[ia],
                                   optimize=paths[ia])
                    out[:] += aux.reshape(out.shape)

            def eval_einsum1(out, eshape, expressions, operands, paths):
                out.reshape(-1)[:] = contract(
                    expressions[0], *operands[0], optimize=paths[0],
                ).reshape(-1)
                for ia in range(1, n_add):
                    out.reshape(-1)[:] += contract(
                        expressions[ia], *operands[ia], optimize=paths[ia],
                    ).reshape(-1)

            def eval_einsum2(out, eshape, expressions, operands, paths):
                out.flat = contract(
                    expressions[0], *operands[0], optimize=paths[0],
                )
                for ia in range(1, n_add):
                    out.ravel()[...] += contract(
                        expressions[ia], *operands[ia], optimize=paths[ia],
                    ).ravel()

            def eval_einsum3(out, eshape, expressions, operands, paths):
                out.ravel()[:] = contract(
                    expressions[0], *operands[0], optimize=paths[0],
                ).ravel()
                for ia in range(1, n_add):
                    out.ravel()[:] += contract(
                        expressions[ia], *operands[ia], optimize=paths[ia],
                    ).ravel()

            def eval_einsum4(out, eshape, expressions, operands, paths):
                vout = out.reshape(eshape)
                contract(expressions[0], *operands[0], out=vout,
                         optimize=paths[0])
                for ia in range(1, n_add):
                    aux = contract(expressions[ia], *operands[ia],
                                   optimize=paths[ia])
                    out[:] += aux.reshape(out.shape)

            eval_fun = self.backend_kwargs.get('eval_fun', 'eval_einsum0')
            eval_einsum = locals()[eval_fun]

        elif self.backend in ('numpy_loop', 'opt_einsum_loop'):
            contract = nm.einsum if self.backend == 'numpy_loop' else oe.contract
            def eval_einsum(out, eshape, expressions, all_slice_ops, paths):
                n_cell = out.shape[0]
                vout = out.reshape(eshape)
                slice_ops = all_slice_ops[0]
                if vout.ndim > 1:
                    for ic in range(n_cell):
                        ops = slice_ops(ic)
                        contract(expressions[0], *ops, out=vout[ic],
                                 optimize=paths[0])

                else: # vout[ic] can be scalar in eval mode.
                    for ic in range(n_cell):
                        ops = slice_ops(ic)
                        vout[ic] = contract(expressions[0], *ops,
                                            optimize=paths[0])

                for ia in range(1, n_add):
                    slice_ops = all_slice_ops[ia]
                    for ic in range(n_cell):
                        ops = slice_ops(ic)
                        vout[ic] += contract(expressions[ia], *ops,
                                             optimize=paths[ia])

        elif self.backend in ('numpy_qloop', 'opt_einsum_qloop'):
            contract = (nm.einsum if self.backend == 'numpy_qloop'
                        else oe.contract)
            def eval_einsum(out, eshape, expressions, all_slice_ops,
                            loop_sizes, paths):
                n_qp = loop_sizes[0]
                vout = out.reshape(eshape)
                slice_ops = all_slice_ops[0]

                ops = slice_ops(0)
                vout[:] = contract(expressions[0], *ops,
                                   optimize=paths[0])
                for iq in range(1, n_qp):
                    ops = slice_ops(iq)
                    vout[:] += contract(expressions[0], *ops,
                                        optimize=paths[0])

                for ia in range(1, n_add):
                    n_qp = loop_sizes[ia]
                    slice_ops = all_slice_ops[ia]
                    for iq in range(n_qp):
                        ops = slice_ops(iq)
                        vout[:] += contract(expressions[ia], *ops,
                                            optimize=paths[ia])

        elif self.backend == 'jax':
            @partial(jax.jit, static_argnums=(0, 1, 2))
            def _eval_einsum(expressions, paths, n_add, operands):
                val = jnp.einsum(expressions[0], *operands[0],
                                 optimize=paths[0])
                for ia in range(1, n_add):
                    val += jnp.einsum(expressions[ia], *operands[ia],
                                      optimize=paths[ia])
                return val

            def eval_einsum(out, eshape, expressions, operands, paths):
                aux = _eval_einsum(expressions, paths, n_add, operands)
                out[:] = nm.asarray(aux.reshape(out.shape))

        elif self.backend == 'jax_vmap':
            def _eval_einsum_cell(expressions, paths, n_add, operands):
                val = jnp.einsum(expressions[0], *operands[0],
                                 optimize=paths[0])
                for ia in range(1, n_add):
                    val += jnp.einsum(expressions[ia], *operands[ia],
                                      optimize=paths[ia])
                return val

            def eval_einsum(out, vmap_eval_cell, eshape, expressions, operands,
                            paths):
                aux = vmap_eval_cell(expressions, paths, n_add,
                                     operands)
                out[:] = nm.asarray(aux.reshape(out.shape))

            eval_einsum = (eval_einsum, _eval_einsum_cell)

        elif self.backend.startswith('dask'):
            scheduler = {'dask_single' : 'single-threaded',
                         'dask_threads' : 'threads'}[self.backend]
            def eval_einsum(out, eshape, expressions, operands, paths):
                _out = da.einsum(expressions[0], *operands[0],
                                 optimize=paths[0])
                for ia in range(1, n_add):
                    aux = da.einsum(expressions[ia],
                                    *operands[ia],
                                    optimize=paths[ia])
                    _out += aux

                out[:] = _out.compute(scheduler=scheduler).reshape(out.shape)

        elif self.backend.startswith('opt_einsum_dask'):
            scheduler = {'opt_einsum_dask_single' : 'single-threaded',
                         'opt_einsum_dask_threads' : 'threads'}[self.backend]

            def eval_einsum(out, eshape, expressions, operands, paths):
                _out = oe.contract(expressions[0], *operands[0],
                                   optimize=paths[0],
                                   backend='dask')
                for ia in range(1, n_add):
                    aux = oe.contract(expressions[ia],
                                      *operands[ia],
                                      optimize=paths[ia],
                                      backend='dask')
                    _out += aux

                out[:] = _out.compute(scheduler=scheduler).reshape(out.shape)

        else:
            raise ValueError('unsupported backend! ({})'.format(self.backend))

        einfo.eval_einsum = eval_einsum

        if self.verbosity:
            output('einsum setup: {} s'.format(timer.stop()))

        return eval_einsum

    def get_operands(self, diff_var):
        einfo = self.einfos[diff_var]
        return get_einsum_ops(einfo.eargs, einfo.ebuilder, self.expr_cache)

    def get_paths(self, expressions, operands):
        memory_limit = self.backend_kwargs.get('memory_limit')

        if ('numpy' in self.backend) or self.backend.startswith('dask'):
            optimize = (self.optimize if memory_limit is None
                        else (self.optimize, memory_limit))
            paths, path_infos = zip(*[nm.einsum_path(
                expressions[ia], *operands[ia],
                optimize=optimize,
            ) for ia in range(len(operands))])

        elif 'opt_einsum' in self.backend:
            paths, path_infos = zip(*[oe.contract_path(
                expressions[ia], *operands[ia],
                optimize=self.optimize,
                memory_limit=memory_limit,
            ) for ia in range(len(operands))])

        elif 'jax' in self.backend:
            paths, path_infos = [], []
            for ia in range(len(operands)):
                path, info = jnp.einsum_path(
                    expressions[ia], *operands[ia],
                    optimize=self.optimize,
                )
                paths.append(tuple(path))
                path_infos.append(info)
            paths = tuple(paths)
            path_infos = tuple(path_infos)

        else:
            raise ValueError('unsupported backend! ({})'.format(self.backend))

        return paths, path_infos

    def get_fargs(self, *args, **kwargs):
        mode, term_mode, diff_var = args[-3:]

        eval_einsum = self.get_function(*args, **kwargs)
        operands = self.get_operands(diff_var)

        einfo = self.einfos[diff_var]
        ebuilder = einfo.ebuilder
        eshape = get_output_shape(ebuilder.out_subscripts[0],
                                  ebuilder.subscripts[0], operands[0])

        out = [eval_einsum, eshape]

        subscripts, operands = ebuilder.apply_layout(
            self.layout, operands, verbosity=self.verbosity,
        )

        self.parsed_expressions = ebuilder.get_expressions(subscripts)
        if self.verbosity:
            output('parsed expressions:', self.parsed_expressions)

        cloop = self.backend in ('numpy_loop', 'opt_einsum_loop', 'jax_vmap')
        qloop = self.backend in ('numpy_qloop', 'opt_einsum_qloop')
        if cloop or qloop:
            loop_index = 'c' if cloop else 'q'
            transform = ebuilder.transform(subscripts, operands,
                                           transformation='loop',
                                           loop_index=loop_index)
            expressions, poperands, all_slice_ops, loop_sizes = transform

            if self.backend == 'jax_vmap':
                all_ics = [get_loop_indices(subs, loop_index)
                           for subs in subscripts]
                vms = (None, None, None, all_ics)
                vmap_eval_cell = jax.jit(jax.vmap(eval_einsum[1], vms, 0),
                                         static_argnums=(0, 1, 2))
                out += [expressions, operands]
                out[:1] = [eval_einsum[0], vmap_eval_cell]

            else:
                out += [expressions, all_slice_ops]
                if qloop:
                    out.append(loop_sizes)

        elif (self.backend.startswith('dask')
              or self.backend.startswith('opt_einsum_dask')):
            c_chunk_size = self.backend_kwargs.get('c_chunk_size')
            da_operands = ebuilder.transform(subscripts, operands,
                                             transformation='dask',
                                             c_chunk_size=c_chunk_size)
            poperands = operands
            expressions = self.parsed_expressions
            out += [expressions, da_operands]

        else:
            poperands = operands
            expressions = self.parsed_expressions
            out += [expressions, operands]

        if einfo.paths is None:
            if self.verbosity > 1:
                ebuilder.print_shapes(subscripts, operands)

            einfo.paths, einfo.path_infos = self.get_paths(
                expressions,
                poperands,
            )
            if self.verbosity > 2:
                for path, path_info in zip(einfo.paths, einfo.path_infos):
                    output('path:', path)
                    output(path_info)

        out += [einfo.paths]

        return out

    def get_eval_shape(self, *args, **kwargs):
        mode, term_mode, diff_var = args[-3:]
        if diff_var is not None:
            raise ValueError('cannot differentiate in {} mode!'
                             .format(mode))

        self.get_function(*args, **kwargs)

        operands = self.get_operands(diff_var)
        ebuilder = self.einfos[diff_var].ebuilder
        out_shape = get_output_shape(ebuilder.out_subscripts[0],
                                     ebuilder.subscripts[0], operands[0])

        dtype = nm.result_type(*operands[0])

        return out_shape, dtype

    def get_normals(self, arg):
        normals = self.get_mapping(arg)[0].normal
        if normals is not None:
            normals = ExpressionArg(name='n({})'.format(arg.name),
                                    arg=normals[..., 0],
                                    kind='ndarray')
        return normals

    def print_expressions(self, diff_var=None):
        operands = self.get_operands(diff_var)

        einfo = self.einfos[diff_var]
        ebuilder = einfo.ebuilder

        subscripts, operands = ebuilder.apply_layout(
            self.layout, operands, verbosity=self.verbosity,
        )

        parsed_expressions = ebuilder.get_expressions(subscripts)
        output('parsed expressions:', parsed_expressions)
        ebuilder.print_shapes(subscripts, operands)

        return ebuilder

class EIntegrateOperatorTerm(ETermBase):
    r"""
    Volume and surface integral of a test function weighted by a scalar
    function :math:`c`.

    :Definition:

    .. math::
        \int_{\cal{D}} q \mbox{ or } \int_{\cal{D}} c q

    :Arguments:
        - material : :math:`c` (optional)
        - virtual  : :math:`q`
    """
    name = 'de_integrate'
    arg_types = ('opt_material', 'virtual')
    arg_shapes = [{'opt_material' : '1, 1', 'virtual' : (1, None)},
                  {'opt_material' : None}]
    integration = ('cell', 'facet')

    def get_function(self, mat, virtual, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        if mat is None:
            fun = self.make_function(
                '0', virtual, mode=mode, diff_var=diff_var,
            )

        else:
            fun = self.make_function(
                '00,0', mat, virtual, mode=mode, diff_var=diff_var,
            )

        return fun

class ELaplaceTerm(ETermBase):
    r"""
    Laplace term with :math:`c` coefficient. Can be
    evaluated. Can use derivatives.

    :Definition:

    .. math::
        \int_{\Omega} \nabla q \cdot \nabla p \mbox{ , }
        \int_{\Omega} c \nabla q \cdot \nabla p

    :Arguments:
        - material: :math:`c`
        - virtual/parameter_1: :math:`q`
        - state/parameter_2: :math:`p`
    """
    name = 'de_laplace'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'parameter_1', 'parameter_2'))
    arg_shapes = [{'opt_material' : '1, 1', 'virtual' : (1, 'state'),
                   'state' : 1, 'parameter_1' : 1, 'parameter_2' : 1},
                  {'opt_material' : None}]
    modes = ('weak', 'eval')

    def get_function(self, mat, virtual, state, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        if mat is None:
            fun = self.make_function(
                '0.j,0.j', virtual, state, mode=mode, diff_var=diff_var,
            )

        else:
            fun = self.make_function(
                '00,0.j,0.j', mat, virtual, state, mode=mode, diff_var=diff_var,
            )

        return fun

class EDotTerm(ETermBase):
    r"""
    Volume and surface :math:`L^2(\Omega)` weighted dot product for both
    scalar and vector fields. Can be evaluated. Can use derivatives.

    :Definition:

    .. math::
        \int_{\cal{D}} q p \mbox{ , } \int_{\cal{D}} \ul{v} \cdot \ul{u}\\
        \int_{\cal{D}} c q p \mbox{ , } \int_{\cal{D}} c \ul{v} \cdot \ul{u}\\
        \int_{\cal{D}} \ul{v} \cdot (\ull{c}\, \ul{u})

    :Arguments:
        - material: :math:`c` or :math:`\ull{c}` (optional)
        - virtual/parameter_1: :math:`q` or :math:`\ul{v}`
        - state/parameter_2: :math:`p` or :math:`\ul{u}`
    """
    name = 'de_dot'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'parameter_1', 'parameter_2'))
    arg_shapes = [{'opt_material' : '1, 1', 'virtual' : (1, 'state'),
                   'state' : 1, 'parameter_1' : 1, 'parameter_2' : 1},
                  {'opt_material' : None},
                  {'opt_material' : '1, 1', 'virtual' : ('D', 'state'),
                   'state' : 'D', 'parameter_1' : 'D', 'parameter_2' : 'D'},
                  {'opt_material' : 'D, D'},
                  {'opt_material' : None},
                  {'opt_material' : '1, 1', 'virtual' : ('N', 'state'),
                   'state' : 'N', 'parameter_1' : 'N', 'parameter_2' : 'N'},
                  {'opt_material' : 'N, N'},
                  {'opt_material' : None}]
    modes = ('weak', 'eval')
    integration = ('cell', 'facet')

    def get_function(self, mat, virtual, state, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        if mat is None:
            fun = self.make_function(
                'i,i', virtual, state, mode=mode, diff_var=diff_var,
            )

        else:
            if mat.shape[-1] > 1:
                fun = self.make_function(
                    'ij,i,j', mat, virtual, state, mode=mode, diff_var=diff_var,
                )

            else:
                fun = self.make_function(
                    '00,i,i', mat, virtual, state, mode=mode, diff_var=diff_var,
                )

        return fun


class EScalarDotMGradScalarTerm(ETermBase):
    r"""
    Volume dot product of a scalar gradient dotted with a material vector with
    a scalar.

    :Definition:

    .. math::
        \int_{\Omega} q \ul{y} \cdot \nabla p \mbox{ , }
        \int_{\Omega} p \ul{y} \cdot \nabla q

    :Arguments 1:
        - material : :math:`\ul{y}`
        - virtual  : :math:`q`
        - state    : :math:`p`

    :Arguments 2:
        - material : :math:`\ul{y}`
        - state    : :math:`p`
        - virtual  : :math:`q`
    """
    name = 'de_s_dot_mgrad_s'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'),
                 ('material', 'parameter_1', 'parameter_2'))
    arg_shapes = [{'material' : 'D, 1',
                   'virtual/grad_state' : (1, None),
                   'state/grad_state' : 1,
                   'virtual/grad_virtual' : (1, None),
                   'state/grad_virtual' : 1,
                   'parameter_1': 1, 'parameter_2': 1}]

    modes = ('grad_state', 'grad_virtual', 'eval')

    def get_function(self, mat, var1, var2, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        return self.make_function(
            'i0,0,0.i', mat, var1, var2, mode=mode, diff_var=diff_var,
        )

class ENonPenetrationPenaltyTerm(ETermBase):
    r"""
    Non-penetration condition in the weak sense using a penalty.

    :Definition:

    .. math::
        \int_{\Gamma} c (\ul{n} \cdot \ul{v}) (\ul{n} \cdot \ul{u})

    :Arguments:
        - material : :math:`c`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """
    name = 'de_non_penetration_p'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'parameter_1', 'parameter_2'))
    arg_shapes = {'material' : '1, 1',
                  'virtual' : ('D', 'state'), 'state' : 'D',
                  'parameter_1' : 'D', 'parameter_2' : 'D'}
    integration = 'facet'
    modes = ('weak', 'eval')

    def get_function(self, mat, virtual, state, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        normals = self.get_normals(state)
        return self.make_function(
            '00,i,i,j,j',
            mat, virtual, normals, state, normals, mode=mode, diff_var=diff_var,
        )

class EDivGradTerm(ETermBase):
    r"""
    Vector field diffusion term.

    :Definition:

    .. math::
        \int_{\Omega} \nabla \ul{v} : \nabla \ul{u} \mbox{ , }
        \int_{\Omega} \nu\ \nabla \ul{v} : \nabla \ul{u}

    :Arguments:
        - material: :math:`\nu` (viscosity, optional)
        - virtual/parameter_1: :math:`\ul{v}`
        - state/parameter_2: :math:`\ul{u}`
    """
    name = 'de_div_grad'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'parameter_1', 'parameter_2'))
    arg_shapes = [{'opt_material' : '1, 1', 'virtual' : ('D', 'state'),
                   'state' : 'D', 'parameter_1' : 'D', 'parameter_2' : 'D'},
                  {'opt_material' : None}]
    modes = ('weak', 'eval')

    def get_function(self, mat, virtual, state, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        if mat is None:
            fun = self.make_function(
                'i.j,i.j', virtual, state, mode=mode, diff_var=diff_var,
            )

        else:
            fun = self.make_function(
                '00,i.j,i.j', mat, virtual, state, mode=mode, diff_var=diff_var,
            )

        return fun


class StokesTractionTerm(ETermBase):
    r"""
    Surface traction term for Stokes flow problem.

    :Definition:

    .. math::
        \int_{\Gamma} \nu \ul{v}\cdot(\nabla \ul{u} \cdot \ul{n})

    :Arguments:
        - material: :math:`\nu` (viscosity, optional)
        - virtual/parameter_1: :math:`\ul{v}`
        - state/parameter_2: :math:`\ul{u}`
    """
    name = 'de_stokes_traction'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'parameter_1', 'parameter_2'))
    arg_shapes = [{'opt_material' : '1, 1', 'virtual' : ('D', 'state'),
                   'state' : 'D', 'parameter_1' : 'D', 'parameter_2' : 'D'},
                  {'opt_material' : None}]
    modes = ('weak', 'eval')
    integration = 'facet_extra'

    def get_function(self, mat, virtual, state, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        normals = self.get_normals(state)
        if mat is None:
            fun = self.make_function(
                'i,i.j,j', virtual, state, normals, mode=mode, diff_var=diff_var,
            )

        else:
            fun = self.make_function(
                '00,i,i.j,j', mat, virtual, state, normals, mode=mode, diff_var=diff_var,
            )

        return fun

class EConvectTerm(ETermBase):
    r"""
    Nonlinear convective term.

    :Definition:

    .. math::
        \int_{\Omega} c ((\ul{u} \cdot \nabla) \ul{u}) \cdot \ul{v}

    :Arguments:
        - material: :math:`c` (optional)
        - virtual/parameter_1: :math:`\ul{v}`
        - state/parameter_2: :math:`\ul{u}`
    """
    name = 'de_convect'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'parameter_1', 'parameter_2'))
    arg_shapes = [{'opt_material' : '1, 1',
                   'virtual' : ('D', 'state'), 'state' : 'D',
                   'parameter_1' : 'D', 'parameter_2' : 'D'},
                  {'opt_material' : None}]
    modes = ('weak', 'eval')

    def get_function(self, coef, virtual, state, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        if coef is None:
            return self.make_function(
                'i,i.j,j', virtual, state, state,
                mode=mode, diff_var=diff_var,
            )

        else:
            return self.make_function(
                '00,i,i.j,j', coef, virtual, state, state,
                mode=mode, diff_var=diff_var,
            )

class ELinearConvectTerm(ETermBase):
    r"""
    Linearized convective term.

    :Definition:

    .. math::
        \int_{\Omega} ((\ul{w} \cdot \nabla) \ul{u}) \cdot \ul{v}

    :Arguments:
        - virtual/parameter_1: :math:`\ul{v}`
        - parameter/parameter_2: :math:`\ul{w}`
        - state/parameter_3: :math:`\ul{u}`
    """
    name = 'de_lin_convect'
    arg_types = (('virtual', 'parameter', 'state'),
                 ('parameter_1', 'parameter_2', 'parameter_3'))
    arg_shapes = {'virtual' : ('D', 'state'), 'parameter' : 'D', 'state' : 'D',
                  'parameter_1' : 'D', 'parameter_2' : 'D', 'parameter_3' : 'D'}
    modes = ('weak', 'eval')

    def get_function(self, virtual, parameter, state, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        return self.make_function(
            'i,i.j,j', virtual, state, parameter, mode=mode, diff_var=diff_var,
        )

class EDivTerm(ETermBase):
    r"""
    Weighted divergence term.

    :Definition:

    .. math::
        \int_{\Omega} \nabla \cdot \ul{v} \mbox { , }
        \int_{\Omega} c \nabla \cdot \ul{v}

    :Arguments:
        - material: :math:`c` (optional)
        - virtual/parameter: :math:`\ul{v}`
    """
    name = 'de_div'
    arg_types = (('opt_material', 'virtual'),
                 ('opt_material', 'parameter'),)
    arg_shapes = [{'opt_material' : '1, 1', 'virtual' : ('D', None),
                   'parameter' : 'D'},
                  {'opt_material' : None}]
    modes = ('weak', 'eval')

    def get_function(self, mat, virtual, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        if mat is None:
            fun = self.make_function(
                'i.i', virtual, mode=mode, diff_var=diff_var,
            )

        else:
            fun = self.make_function(
                '00,i.i', mat, virtual, mode=mode, diff_var=diff_var,
            )

        return fun


class EGradTerm(ETermBase):
    r"""
    Weighted gradient term.

    :Definition:

    .. math::
        \int_{\Omega} \nabla \ul{v} \mbox { , }
        \int_{\Omega} c \nabla \ul{v} \mbox { , }
        \int_{\Omega} \ul{c} \cdot \nabla \ul{v} \mbox { , }
        \int_{\Omega} \ull{c} \cdot \nabla \ul{v}

    :Arguments:
        - material: :math:`c` (optional)
        - virtual/parameter: :math:`\ul{v}`
    """
    name = 'de_grad'
    arg_types = ('opt_material', 'parameter')
    arg_shapes = [{'opt_material' : '1, 1', 'parameter' : 'N'},
                  {'opt_material' : 'N, 1'}, {'opt_material' : 'N, N'},
                  {'opt_material' : None}]

    def get_function(self, mat, virtual, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        if mat is None:
            fun = self.make_function(
                'i.j', virtual, mode=mode, diff_var=diff_var,
            )

        else:
            dim = virtual.dim
            if mat.shape[-2:] == (1, 1):
                expr = '00,i.j'
            elif mat.shape[-2:] in [(dim, 1), (dim, dim)]:
                expr = 'jk,i.j'
            else:
                raise ValueError(f'wrong matrial shape! ({mat.shape[-2:]})')

            fun = self.make_function(
                expr, mat, virtual, mode=mode, diff_var=diff_var,
            )

        return fun


class EStokesTerm(ETermBase):
    r"""
    Stokes problem coupling term. Corresponds to weak forms of gradient and
    divergence terms.

    :Definition:

    .. math::
        \int_{\Omega} p\, \nabla \cdot \ul{v} \mbox{ , }
        \int_{\Omega} q\, \nabla \cdot \ul{u}\\
        \int_{\Omega} c\, p\, \nabla \cdot \ul{v} \mbox{ , }
        \int_{\Omega} c\, q\, \nabla \cdot \ul{u}

    :Arguments 1:
        - material: :math:`c` (optional)
        - virtual/parameter_v: :math:`\ul{v}`
        - state/parameter_s: :math:`p`

    :Arguments 2:
        - material : :math:`c` (optional)
        - state    : :math:`\ul{u}`
        - virtual  : :math:`q`
    """
    name = 'de_stokes'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'state', 'virtual'),
                 ('opt_material', 'parameter_v', 'parameter_s'))
    arg_shapes = [{'opt_material' : '1, 1',
                   'virtual/grad' : ('D', None), 'state/grad' : 1,
                   'virtual/div' : (1, None), 'state/div' : 'D',
                   'parameter_v' : 'D', 'parameter_s' : 1},
                  {'opt_material' : None}]
    modes = ('grad', 'div', 'eval')

    def get_function(self, coef, vvar, svar, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        if coef is None:
            fun = self.make_function(
                'i.i,0', vvar, svar, mode=mode, diff_var=diff_var,
            )

        else:
            fun = self.make_function(
                '00,i.i,0', coef, vvar, svar, mode=mode, diff_var=diff_var,
            )

        return fun

class ELinearElasticTerm(ETermBase):
    r"""
    General linear elasticity term, with :math:`D_{ijkl}` given in
    the usual matrix form exploiting symmetry: in 3D it is :math:`6\times6`
    with the indices ordered as :math:`[11, 22, 33, 12, 13, 23]`, in 2D it is
    :math:`3\times3` with the indices ordered as :math:`[11, 22, 12]`.

    :Definition:

    .. math::
        \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})

    :Arguments:
        - material: :math:`D_{ijkl}`
        - virtual/parameter_1: :math:`\ul{v}`
        - state/parameter_2: :math:`\ul{u}`
    """
    name = 'de_lin_elastic'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'parameter_1', 'parameter_2'))
    arg_shapes = {'material' : 'S, S', 'virtual' : ('D', 'state'),
                  'state' : 'D', 'parameter_1' : 'D', 'parameter_2' : 'D'}
    modes = ('weak', 'eval')

    def get_function(self, mat, virtual, state, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        return self.make_function(
            'IK,s(i:j)->I,s(k:l)->K', mat, virtual, state,
            mode=mode, diff_var=diff_var,
        )

class ECauchyStressTerm(ETermBase):
    r"""
    Evaluate Cauchy stress tensor.

    It is given in the usual vector form exploiting symmetry: in 3D it has 6
    components with the indices ordered as :math:`[11, 22, 33, 12, 13, 23]`, in
    2D it has 3 components with the indices ordered as :math:`[11, 22, 12]`.

    :Definition:

    .. math::
        \int_{\Omega} D_{ijkl} e_{kl}(\ul{w})

    :Arguments:
        - material  : :math:`D_{ijkl}`
        - parameter : :math:`\ul{w}`
    """
    name = 'de_cauchy_stress'
    arg_types = ('material', 'parameter')
    arg_shapes = {'material' : 'S, S', 'parameter' : 'D'}

    def get_function(self, mat, parameter, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        return self.make_function(
            'IK,s(k:l)->K', mat, parameter, mode=mode, diff_var=diff_var,
        )


class ENonSymElasticTerm(ETermBase):
    r"""
    Elasticity term with non-symmetric gradient. The indices of matrix
    :math:`D_{ijkl}` are ordered as
    :math:`[11, 12, 13, 21, 22, 23, 31, 32, 33]` in 3D and as
    :math:`[11, 12, 21, 22]` in 2D.

    :Definition:

    .. math::
        \int_{\Omega} \ull{D} \nabla \ul{v} : \nabla \ul{u}

    :Arguments:
        - material: :math:`\ull{D}`
        - virtual/parameter_1: :math:`\ul{v}`
        - state/parameter_2: :math:`\ul{u}`
    """
    name = 'de_nonsym_elastic'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'parameter_1', 'parameter_2'))
    arg_shapes = {'material': 'D2, D2', 'virtual': ('D', 'state'),
                  'state': 'D', 'parameter_1': 'D', 'parameter_2': 'D'}
    modes = ('weak', 'eval')

    def get_function(self, mat, virtual, state, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        fun = self.make_function(
            'IK,v(i.j)->I,v(k.l)->K', mat, virtual, state,
            mode=mode, diff_var=diff_var,
        )

        return fun


def sym2nonsym(sym_obj, axes=[3]):
    nonsym_tab = {
        3: nm.array([0, 2, 2, 1]),
        6: nm.array([0, 3, 4, 3, 1, 5, 4, 5, 2]),
    }

    if isinstance(axes, int):
        axes = [axes]

    sh = sym_obj.shape
    sym = sh[axes[0]]
    idxs = nonsym_tab[sym]
    ix = (range(sh[k]) if k not in axes else idxs for k in range(len(sh)))

    return sym_obj[nm.ix_(*ix)]


class EDiffusionTerm(ETermBase):
    r"""
    General diffusion term.

    :Definition:

    .. math::
        \int_{\Omega} K_{ij} \nabla_i q\, \nabla_j p

    :Arguments:
        - material: :math:`K_{ij}`
        - virtual/parameter_1: :math:`q`
        - state/parameter_2: :math:`p`
    """
    name = 'de_diffusion'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'parameter_1', 'parameter_2'))
    arg_shapes = {'material' : 'D, D', 'virtual': (1, 'state'), 'state': 1,
                  'parameter_1': 1, 'parameter_2': 1}
    modes = ('weak', 'eval')

    def get_function(self, mat, vvar, svar,
                     mode=None, term_mode=None, diff_var=None, **kwargs):

        fun = self.make_function(
            'ij,0.i,0.j', mat, vvar, svar, mode=mode, diff_var=diff_var
        )

        return fun


class ELinearTractionTerm(ETermBase):
    r"""
    Linear traction term. The material parameter can have one of the
    following shapes:

        - 1 or (1, 1) - a given scalar pressure
        - (D, 1) - a traction vector
        - (S, 1) or (D, D) - a given stress in symmetric or
          non-symmetric tensor storage (in symmetric storage indicies are order
          as follows: 2D: [11, 22, 12], 3D: [11, 22, 33, 12, 13, 23])

    :Definition:

    .. math::
        \int_{\Gamma} \ul{v} \cdot \ul{n}
            \mbox{ , }
        \int_{\Gamma} c\, \ul{v} \cdot \ul{n}\\
        \int_{\Gamma} \ul{v} \cdot (\ull{\sigma}\, \ul{n})
            \mbox{ , }
        \int_{\Gamma} \ul{v} \cdot \ul{f}

    :Arguments:
        - material: :math:`c`, :math:`\ul{f}`, :math:`\ul{\sigma}`
          or :math:`\ull{\sigma}`
        - virtual/parameter: :math:`\ul{v}`
    """
    name = 'de_surface_ltr'
    arg_types = (('opt_material', 'virtual'),
                 ('opt_material', 'parameter'))
    arg_shapes = [{'opt_material': 'S, 1', 'virtual': ('D', None),
                   'parameter': 'D'},
                  {'opt_material': None}, {'opt_material': '1, 1'},
                  {'opt_material': 'D, 1'}, {'opt_material': 'D, D'}]
    modes = ('weak', 'eval')
    integration = 'facet'

    def get_function(self, traction, vvar, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        sg, _ = self.get_mapping(vvar)
        nv = sg.normal
        _, n_qp, dim, _ = nv.shape

        tdim, tdim2 = (None, None) if traction is None else traction.shape[2:]

        if tdim == dim and tdim2 == 1:  # force vector
            force = traction
        else:
            sym = (dim + 1) * dim // 2

            if tdim is None:
                force = nv
            elif tdim == 1:  # scalar value
                force = traction * nv
            elif tdim == dim and tdim2 == dim:  # traction tensor - full
                force = dot_sequences(traction, nv, 'AB')
            elif tdim == dim**2 and tdim2 == 1:  # traction tensor - full, vector
                aux = traction.reshape((-1, n_qp, dim, dim))
                force = dot_sequences(aux, nv, 'AB')
            elif tdim == sym and tdim2 == 1:  # traction tensor - symetric
                aux = sym2nonsym(traction, [2]).reshape((-1, n_qp, dim, dim))
                force = dot_sequences(aux, nv, 'AB')

        fun = self.make_function(
            'i,i', (force[..., 0], 'opt_material'), vvar,
            mode=mode, diff_var=diff_var,
        )

        return fun

class SurfaceFluxOperatorTerm(ETermBase):
    r"""
    Surface flux operator term.

    :Definition:

    .. math::
        \int_{\Gamma} q \ul{n} \cdot \ull{K} \cdot \nabla p \mbox{ , }
        \int_{\Gamma} p \ul{n} \cdot \ull{K} \cdot \nabla q

    :Arguments 1:
        - material : :math:`\ull{K}`
        - virtual  : :math:`q`
        - state    : :math:`p`

    :Arguments 2:
        - material : :math:`\ull{K}`
        - state    : :math:`p`
        - virtual  : :math:`q`
    """
    name = 'de_surface_flux'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'),
                 ('material', 'parameter_1', 'parameter_2'))
    arg_shapes = [{'material' : 'D, D',
                   'virtual/grad_state' : (1, None),
                   'state/grad_state' : 1,
                   'virtual/grad_virtual' : (1, None),
                   'state/grad_virtual' : 1,
                   'parameter_1': 1, 'parameter_2': 1}]
    integration = 'facet_extra'
    modes = ('grad_state', 'grad_virtual', 'eval')

    def get_function(self, mat, var1, var2, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        normals = self.get_normals(var2)
        return self.make_function(
            'i,ij,0,0.j',
            normals, mat, var1, var2, mode=mode, diff_var=diff_var,
        )

class SurfacePiezoFluxOperatorTerm(ETermBase):
    r"""
    Surface piezoelectric flux operator term.

    Corresponds to the electric flux due to mechanically induced electrical
    displacement.

    :Definition:

    .. math::
        \int_{\Gamma} q g_{kij} e_{ij}(\ul{u}) n_k \mbox{ , }
        \int_{\Gamma} p g_{kij} e_{ij}(\ul{v}) n_k

    :Arguments 1:
        - material            : :math:`g_{kij}`
        - virtual/parameter_1 : :math:`q`
        - state/parameter_2   : :math:`\ul{u}`

    :Arguments 2:
        - material : :math:`g_{kij}`
        - state    : :math:`p`
        - virtual  : :math:`\ul{v}`
    """
    name = 'de_surface_piezo_flux'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'),
                 ('material', 'parameter_1', 'parameter_2'))
    arg_shapes = [{'material' : 'D, S',
                   'virtual/grad_state' : (1, None),
                   'state/grad_state' : 'D',
                   'virtual/grad_virtual' : ('D', None),
                   'state/grad_virtual' : 1,
                   'parameter_1': 1, 'parameter_2': 'D'}]
    integration = 'facet_extra'
    modes = ('grad_state', 'grad_virtual', 'eval')

    def get_function(self, mat, var1, var2, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        normals = self.get_normals(var2)
        return self.make_function(
            'i,iJ,0,s(a:b)->J',
            normals, mat, var1, var2, mode=mode, diff_var=diff_var,
        )

class SurfacePiezoFluxTerm(ETermBase):
    r"""
    Surface piezoelectric flux term.

    Corresponds to the electric flux due to mechanically induced electrical
    displacement.

    :Definition:

    .. math::
        \int_{\Gamma} g_{kij} e_{ij}(\ul{u}) n_k

    :Arguments 1:
        - material  : :math:`g_{kij}`
        - parameter : :math:`\ul{u}`
    """
    name = 'ev_surface_piezo_flux'
    arg_types = ('material', 'parameter')
    arg_shapes = {'material' : 'D, S', 'parameter': 'D'}
    integration = 'facet_extra'
    modes = ('eval',)

    def get_function(self, mat, var1, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        normals = self.get_normals(var1)
        return self.make_function(
            'i,iJ,s(a:b)->J',
            normals, mat, var1, mode=mode, diff_var=diff_var,
        )

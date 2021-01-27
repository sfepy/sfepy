from __future__ import absolute_import
from copy import copy

import numpy as nm

from sfepy.base.testing import TestCommon
from sfepy.base.base import ordered_iteritems
from sfepy import data_dir

filename_meshes = [data_dir + '/meshes/elements/%s_2.mesh' % geom
                   for geom in ['1_2', '2_3', '2_4', '3_4', '3_8', '3_2_4']]

def make_term_args(arg_shapes, arg_kinds, arg_types, ats_mode, domain,
                   material_value=None, poly_space_base=None):
    from sfepy.base.base import basestr
    from sfepy.discrete import FieldVariable, Material, Variables, Materials
    from sfepy.discrete.fem import Field
    from sfepy.solvers.ts import TimeStepper
    from sfepy.mechanics.tensors import dim2sym

    omega = domain.regions['Omega']
    dim = domain.shape.dim
    sym = dim2sym(dim)

    def _parse_scalar_shape(sh):
        if isinstance(sh, basestr):
            if sh == 'D':
                return dim

            elif sh == 'D2':
                return dim**2

            elif sh == 'S':
                return sym

            elif sh == 'N': # General number ;)
                return 1

            else:
                return int(sh)

        else:
            return sh

    def _parse_tuple_shape(sh):
        if isinstance(sh, basestr):
            return [_parse_scalar_shape(ii.strip()) for ii in sh.split(',')]

        else:
            return (int(sh),)

    args = {}
    str_args = []
    materials = []
    variables = []
    for ii, arg_kind in enumerate(arg_kinds):
        if arg_kind != 'ts':
            if ats_mode is not None:
                extended_ats = arg_types[ii] + ('/%s' % ats_mode)

            else:
                extended_ats = arg_types[ii]

            try:
                sh = arg_shapes[arg_types[ii]]

            except KeyError:
                sh = arg_shapes[extended_ats]

        if arg_kind.endswith('variable'):
            shape = _parse_scalar_shape(sh[0] if isinstance(sh, tuple) else sh)
            field = Field.from_args('f%d' % ii, nm.float64, shape, omega,
                                    approx_order=1,
                                    poly_space_base=poly_space_base)

            if arg_kind == 'virtual_variable':
                if sh[1] is not None:
                    istate = arg_types.index(sh[1])

                else:
                    # Only virtual variable in arguments.
                    istate = -1
                    # -> Make fake variable.
                    var = FieldVariable('u-1', 'unknown', field)
                    var.set_constant(0.0)
                    variables.append(var)

                var = FieldVariable('v', 'test', field,
                                    primary_var_name='u%d' % istate)

            elif arg_kind == 'state_variable':
                var = FieldVariable('u%d' % ii, 'unknown', field)
                var.set_constant(0.0)

            elif arg_kind == 'parameter_variable':
                var = FieldVariable('p%d' % ii, 'parameter', field,
                                    primary_var_name='(set-to-None)')
                var.set_constant(0.0)

            variables.append(var)
            str_args.append(var.name)
            args[var.name] = var

        elif arg_kind.endswith('material'):
            if sh is None: # Switched-off opt_material.
                continue

            prefix = ''
            if isinstance(sh, basestr):
                aux = sh.split(':')
                if len(aux) == 2:
                    prefix, sh = aux

            if material_value is None:
                material_value = 1.0

            shape = _parse_tuple_shape(sh)
            if (len(shape) > 1) or (shape[0] > 1):
                if ((len(shape) == 2) and (shape[0] ==  shape[1])
                    and (material_value != 0.0)):
                    # Identity matrix.
                    val = nm.eye(shape[0], dtype=nm.float64)

                else:
                    # Array.
                    val = nm.empty(shape, dtype=nm.float64)
                    val.fill(material_value)

                values = {'%sc%d' % (prefix, ii)
                          : val}

            elif (len(shape) == 1) and (shape[0] == 1):
                # Single scalar as a special value.
                values = {'.c%d' % ii : material_value}

            else:
                raise ValueError('wrong material shape! (%s)' % shape)

            mat = Material('m%d' % ii, values=values)

            materials.append(mat)
            str_args.append(mat.name + '.' + 'c%d' % ii)
            args[mat.name] = mat

        elif arg_kind == 'ts':
            ts = TimeStepper(0.0, 1.0, 1.0, 5)
            str_args.append('ts')
            args['ts'] = ts

        else:
            str_args.append('user%d' % ii)
            args[str_args[-1]] = None

    materials = Materials(materials)
    variables = Variables(variables)

    return args, str_args, materials, variables

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        from sfepy.discrete import Integral
        from sfepy.discrete.fem import Mesh, FEDomain

        domains = []
        for filename in filename_meshes:
            mesh = Mesh.from_file(filename)
            domain = FEDomain('domain_%s' % mesh.name.replace(data_dir, ''),
                              mesh)
            domain.create_region('Omega', 'all')
            domain.create_region('Gamma', 'vertices of surface', 'facet')

            domains.append(domain)

        integral = Integral('i', order=3)
        qp_coors, qp_weights = integral.get_qp('3_8')
        custom_integral = Integral('i', coors=qp_coors, weights=qp_weights,
                                   order='custom')

        test = Test(domains=domains, integral=integral,
                    custom_integral=custom_integral,
                    conf=conf, options=options)
        return test

    def test_term_call_modes(self):
        from sfepy.terms import term_table
        ok = True

        failed = []
        for domain in self.domains:
            self.report('domain: %s' % domain.name)

            domain_geometry = list(domain.geom_els.values())[0].name
            if domain.shape.dim != domain.shape.tdim:
                domain_geometry = '%d_%s' % (domain.shape.dim, domain_geometry)

            for _, term_cls in ordered_iteritems(term_table):
                if (domain_geometry not in term_cls.geometries) \
                   or ("dg" in term_cls.name) \
                   or (term_cls.name == "dw_ns_dot_grad_s"):
                    continue

                vint = ('volume', 'point', 'custom')
                rname = 'Omega' if term_cls.integration in vint else 'Gamma'

                self.report('<-- %s ...' % term_cls.name)

                if rname == 'Gamma' and domain.mesh.dim == 1:
                    self.report('--> 1D Gamma region: not tested!')

                elif term_cls.arg_shapes:
                    try:
                        _ok = self._test_single_term(term_cls, domain, rname)

                    except:
                        _ok = False

                    if not _ok:
                        failed.append((domain.name, term_cls.name))

                    ok = ok and _ok
                    self.report('--> ok: %s' % _ok)

                else:
                    self.report('--> not tested!')

        self.report('failed:', failed)

        return ok

    def _test_single_term(self, term_cls, domain, rname):
        from sfepy.terms import Term
        from sfepy.terms.terms import get_arg_kinds

        ok = True

        term_call = term_cls.name + '(%s)'

        arg_shapes_list = term_cls.arg_shapes
        if not isinstance(arg_shapes_list, list):
            arg_shapes_list = [arg_shapes_list]

        if term_cls.integration != 'custom':
            integral = self.integral

        else:
            integral = self.custom_integral

        poly_space_base = getattr(term_cls, 'poly_space_base', 'lagrange')

        prev_shapes = {}
        for _arg_shapes in arg_shapes_list:
            # Unset shapes are taken from the previous iteration.
            arg_shapes = copy(prev_shapes)
            arg_shapes.update(_arg_shapes)
            prev_shapes = arg_shapes

            self.report('arg_shapes:', arg_shapes)
            arg_types = term_cls.arg_types
            if not isinstance(arg_types[0], tuple):
                arg_types = (arg_types,)

            for iat, ats in enumerate(arg_types):
                self.report('arg_types:', ats)

                arg_kinds = get_arg_kinds(ats)
                modes = getattr(term_cls, 'modes', None)
                mode = modes[iat] if modes is not None else None

                if 'dw_s_dot_grad_i_s' in term_cls.name:
                    material_value = 0.0

                else:
                    material_value = 1.0
                aux = make_term_args(arg_shapes, arg_kinds, ats, mode, domain,
                                     material_value=material_value,
                                     poly_space_base=poly_space_base)
                args, str_args, materials, variables = aux

                self.report('args:', str_args)

                name = term_call % (', '.join(str_args))
                term = Term.new(name, integral, domain.regions[rname], **args)
                term.setup()

                call_mode = 'weak' if term.names.virtual else 'eval'
                self.report('call mode:', call_mode)

                out = term.evaluate(mode=call_mode, ret_status=True)

                if call_mode == 'eval':
                    vals, status = out
                    vals = nm.array(vals)

                else:
                    vals, iels, status = out

                if isinstance(vals, tuple):
                    # Dynamic connectivity terms.
                    vals = vals[0]

                _ok = nm.isfinite(vals).all()
                ok = ok and _ok
                self.report('values shape: %s' % (vals.shape,))
                if not _ok:
                    self.report('values are not finite!')
                    self.report(vals)

                _ok = status == 0
                if not _ok:
                    self.report('status is %s!' % status)

                ok = ok and _ok

                if term.names.virtual:
                    # Test differentiation w.r.t. state variables in the weak
                    # mode.
                    svars = term.get_state_variables(unknown_only=True)
                    for svar in svars:
                        vals, iels, status = term.evaluate(mode=call_mode,
                                                           diff_var=svar.name,
                                                           ret_status=True)
                        if isinstance(vals, tuple):
                            # Dynamic connectivity terms.
                            vals = vals[0]

                        _ok = status == 0
                        ok = ok and _ok
                        self.report('diff: %s' % svar.name)
                        if not _ok:
                            self.report('status is %s!' % status)

                        _ok = nm.isfinite(vals).all()
                        ok = ok and _ok
                        self.report('values shape: %s' % (vals.shape,))
                        if not _ok:
                            self.report('values are not finite!')
                            self.report(vals)

        return ok

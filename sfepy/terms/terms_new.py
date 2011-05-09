"""
todo:
    - get row variable, col variable (if diff_var)
    - determine out shape
    - set current group to all variable arguments
    - loop over row/col dofs:
        - call term

(*) for given group, Variable has to be able to:
    - evaluate value in quadrature points
    - evaluate gradient in quadrature points
    - ? evaluate divergence in quadrature points
    - ? lazy evaluation, cache!

?? base function gradients in space elements stored now in terms - in
geometry, shared dict of geometries belongs to Equations
-> where to cache stuff?

"""
import numpy as nm

from sfepy.base.base import output
from sfepy.terms.terms import Term
from sfepy.linalg import dot_sequences as dot

class NewTerm(Term):

    def get_shape_kind(self):
        if self.integration == 'surface':
            shape_kind = 'surface'

        elif self.integration in ('volume', 'surface_extra'):
            shape_kind = 'volume'

        else:
            raise NotImplementedError('unsupported term integration! (%s)'
                                      % self.integration)

        return shape_kind

    def evaluate(self, mode='eval', diff_var=None, **kwargs):
        shape_kind = self.get_shape_kind()

        if mode == 'eval':
            val = 0.0
            for ig in self.iter_groups():
                args = self.get_args(**kwargs)
                for aux, iels, status in self(*args,
                                              call_mode='d_eval', **kwargs):
                    val += self.sign * aux

        elif mode in ('el_avg', 'qp'):
            pass

        elif mode == 'weak':
            varr = self.get_virtual_variable()

            vals = []
            iels = []

            if diff_var is None:
                for ig in self.iter_groups():
                    args = self.get_args(**kwargs)

                    aux = varr.get_data_shape(ig, self.integral,
                                              shape_kind, self.region.name)
                    n_elr, n_qpr, n_cr, n_enr = aux
                    n_row = n_cr * n_enr

                    shape = (n_elr, 1, n_row, 1)
                    val = nm.zeros(shape, dtype=varr.dtype)
                    for ir in varr.iter_dofs():
                        try:
                            val_qp = self(*args, **kwargs)

                        except ValueError:
                            output('%s term evaluation failed!' % self.name)
                            raise

                        _val = self.integrate(val_qp, varr)
                        val[..., ir, 0] = _val


                    val.append(self.sign * val)
                    iels.append((ig, nm.arange(n_elr, dtype=nm.int32)))

            else:
                varc = self.get_variables()[diff_var]

                for ig in self.iter_groups():
                    args = self.get_args(**kwargs)

                    aux = varr.get_data_shape(ig, self.integral,
                                              shape_kind, self.region.name)
                    n_elr, n_qpr, n_cr, n_enr = aux
                    n_row = n_cr * n_enr

                    aux = varc.get_data_shape(ig, self.integral,
                                              shape_kind, self.region.name)
                    n_elc, n_qpc, n_cc, n_enc = aux
                    n_col = n_cr * n_enr

                    shape = (n_elr, 1, n_row, n_col)
                    val = nm.zeros(shape, dtype=varr.dtype)
                    for ir in varr.iter_dofs():
                        for ic in varc.iter_dofs():
                            try:
                                val_qp = self(*args, **kwargs)

                            except ValueError:
                                output('%s term evaluation failed!' % self.name)
                                raise

                        _val = self.integrate(val_qp, varr)
                        val[..., ir, ic] = _val


                    val.append(self.sign * val)
                    iels.append((ig, nm.arange(n_elr, dtype=nm.int32)))

        # Setup return value.
        if mode == 'eval':
            out = (vals,)

        else:
            out = (vals, iels)

        if len(out) == 1:
            out = out[0]

        # Hack: add zero status.
        return out + (0,)

    def set_current_group(self, ig):
        """
        Set current group for the term and all variable in its
        arguments.
        """
        self.char_fun.set_current_group(ig)

        shape_kind = self.get_shape_kind()
        for var in self.get_variables():
            var.set_current_group(ig, self.integral,
                                  shape_kind, self.region.name)

class NewDiffusionTerm(NewTerm):
    """
    """
    name = 'dw_new_diffusion'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'parameter_1', 'parameter_2'))
    modes = ('weak', 'eval')

    def __call__(self, mat, virtual, state, **kwargs):

        val = mat * dot(virtual.grad(), state.grad(), 'ATB')

        return val

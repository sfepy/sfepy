import numpy as np

from sfepy.base.base import (as_float_or_complex, get_default, assert_,
                             Container, Struct, basestr, goptions)

from sfepy.terms.terms import Term
from scipy.sparse import lil_matrix, csr_matrix

class IntFE(Term):
    r"""
    General diffusion term for integration over finite elements.
    """
    name = 'intFE'

    def __init__(self, name, arg_str, integral, region, **kwargs):
        """
        It is required here to pass the other term parameters (fun_args, etc.)
        to self.desc.
        """
        self.name = name
        self.arg_str = arg_str
        self.region = region
        self.desc = kwargs['desc']
        self._integration = self.integration
        self.sign = 1.0

        self.set_integral(integral)

    def setup(self):
        """
        It is here because of self.setup_formal_args() which do not fit to
        arg_methods and because of self._kwargs which is not defined here
        (functionality ???).
        """
        from sfepy.terms.terms import CharacteristicFunction
        self.char_fun = CharacteristicFunction(self.region)
        self.function = Struct.get(self, 'function', None)

        self.step = 0
        self.dt = 1.0
        self.is_quasistatic = False
        self.has_region = True

        if not hasattr(self, 'arg_names'):
            self.setup_formal_args()

    def setup_formal_args(self):
        """
        Auxiliary function that setup arguments of integral.
        """
        self.arg_names = []
        self.arg_steps = {}
        self.arg_traces = {}
        self.arg_methods = []
        arg_list = self.arg_str.split(',')
        for arg in arg_list:
            arg_split = arg.split('.')
            arg_name = arg_split[0].strip()
            self.arg_names.append(arg_name)
            if len(arg_split) < 2:
                self.arg_methods.append(['val'])
            else:
                arg_method = []
                for arg_method_str in arg_split[1:]:
                    arg_method.append(arg_method_str.strip())
                self.arg_methods.append(arg_method)

            self.arg_steps[arg_name] = 0
            self.arg_traces[arg_name] = False

    def assign_args(self, variables, materials, user=None):
        """
        Check term argument existence in variables, materials, user data
        and assign the arguments to terms.
        """
        if user is None:
            user = {}

        self.kwargs = {}
        self.arg_vars = []
        self.arg_vars_method = []
        self.arg_mats = []
        self.names = Struct(name='arg_names',
                            material=[], material_name=[],
                            variable=[], variable_method=[],
                            user=[], state=[], virtual=[], parameter=[])

        for ii, arg_name in enumerate(self.arg_names):

            if isinstance(arg_name, basestr):

                if arg_name in variables.names:
                    self.kwargs[arg_name] = variables[arg_name]
                    self.arg_vars.append(variables[arg_name])
                    self.arg_vars_method.append(self.arg_methods[ii])
                    self.names.variable.append(arg_name)
                    if self.arg_vars[-1].kind == 'unknown':
                        self.names.state.append(arg_name)
                    elif self.arg_vars[-1].kind == 'test':
                        self.names.virtual.append(arg_name)
                    else:
                        msg = 'Variable kind should be unknown or test but \
                            it is %s' % (str(self.arg_vars[-1].kind))
                        NotImplementedError(msg)

                elif arg_name in materials.names:
                    self.arg_methods[ii] = self.arg_methods[ii][0]
                    self.kwargs[arg_name] = materials[arg_name]
                    self.arg_mats.append(materials[arg_name])
                    self.names.material.append((arg_name,
                                                self.arg_methods[ii]))
                    self.names.material_name.append(arg_name)

                elif arg_name in user:
                    self.kwargs[arg_name] = user[arg_name]

                else:
                    raise ValueError('argument %s not found!' % arg_name)

        self.args = []
        for arg_name in self.arg_names:
            self.args.append(self.kwargs[arg_name])

        self.classify_args()
        self.check_args()

    def classify_args(self):
        """
        Classify types of the term arguments and find matching call
        signature.

        A state variable can be in place of a parameter variable and
        vice versa.
        """
        self.n_virtual = len(self.names.virtual)
        if self.n_virtual > 1:
            raise ValueError('at most one virtual variable is allowed! (%d)'
                             % self.n_virtual)

        self.setup_integration()

    def check_args(self):
        pass

    def get_materials(self, join=False):
        materials = self.arg_mats

        if join:
            materials = list(set(materials))

        return materials

    def evaluate(self, asm_obj=None, mode='eval', dw_mode='vector',
                 diff_var=None, standalone=True, ret_status=False, **kwargs):
        """
        This function evaluates the term via self.eval_lin method. The standard
        parameters such as mode, dw_mode are processed and passed to
        self.eval_lin method.

        Parameters
        ----------
        mode : 'eval' (default), or 'weak'
            The term evaluation mode.
        dw_mode : 'vector' or 'matrix'
            It distinguish between evaluation a linear or bilinear form
            leading to vector or matrix assemble object.

        Returns
        -------
        val : float or array
            In 'eval' mode, the term returns a single value (the
            integral, it does not need to be a scalar), while in 'weak'
            mode it returns an array for each element.
        status : int, optional
            The flag indicating evaluation success (0) or failure
            (nonzero). Only provided if `ret_status` is True.
        iels : array of ints, optional
            The local elements indices in 'weak' mode. Only provided in
            non-'eval' modes.
        """

        if mode == 'eval':
            asm_obj = self.eval_lin(asm_obj=asm_obj, active_dof=False)

        elif mode == 'weak':

            vvar = self.get_virtual_variable()

            if dw_mode == 'vector':
                asm_obj = self.eval_lin(asm_obj=asm_obj, linear=[vvar.name],
                                        active_dof=True)

            elif dw_mode == 'matrix':
                svar = self.get_state_variables()
                svar_name = svar[0].name
                asm_obj = self.eval_lin(asm_obj=asm_obj,
                                        linear=[vvar.name, svar_name],
                                        active_dof=True)
            else:
                NotImplementedError('This type of dw_mode (%s) is not \
                                    supported.' % (dw_mode))

        else:
            NotImplementedError('This type of mode (%s) is not \
                                    supported.' % (mode))

        return asm_obj

    def eval_lin(self, asm_obj=None, linear=[], active_dof=False,
                 global_matrix=True):
        r"""
        The key method that evaluates integral over elements and assemble it
        to a predefined object.
        Integral arguments are sorted to: material, parameter, and linear.
        The number of linear variables determines if the integral represents
        functional, linear form, or bilinear form corresponding also
        to assemble object, which is scalar, vector, or matrix respectively.

        Parameters
        ----------
        asm_obj : numpy.array
            the object for assembling to; it can be one of scalar, vector,
            or matrix
        linear : list of string
            variables that are linear in integral;
            no linear variable represent functional returns scalar,
            one linear variable returns vector and represents linear form,
            two linear variables returns matrix and represents bilinear form
        active_dof : boolean
            if True it will return the submatrix according to
            active DOFs

        Returns
        -------
        asm_obj : numpy.array
            assemble object
        """
        kwargs = self.kwargs
        args = self.args
        arg_vars = self.arg_vars
        arg_mats = self.arg_mats
        n_args = len(arg_vars) + len(arg_mats)
        region = self.region

        # sorting of arguments to parameter and linear
        vars_lin = []; vars_par = []; vars_lin_names = []; vars_par_names = []
        for arg in self.arg_vars:
            if arg.name in linear: # linear argument
                vars_lin.append(arg)
                vars_lin_names.append(arg.name)
            else: # parameter argument
                vars_par.append(arg)
                vars_par_names.append(arg.name)

        # create assemble object for term
        n_args_lin = len(vars_lin_names)
        if n_args_lin == 0: # scalar
            asm_obj_T = np.zeros(1)
            val_loc_shape = 1

            def asm_fun(asm_obj_T, econn_el, val):
                asm_obj_T[0] += val

        elif n_args_lin == 1: # vector
            asm_obj_T = np.zeros(kwargs[vars_lin_names[0]].n_dof)

            def asm_fun(asm_obj_T, econn_el, val):
                asm_obj_T[econn_el[0][ilin[0]]] += val

        elif n_args_lin == 2: # matrix
            shape = (kwargs[vars_lin_names[0]].n_dof,
                     kwargs[vars_lin_names[1]].n_dof)
            asm_obj_T = lil_matrix(shape)

            def asm_fun(asm_obj_T, econn_el, val):
                asm_obj_T[econn_el[0][ilin[0]], econn_el[1][ilin[1]]] += val

        coors = region.domain.mesh.coors # coordinates of nodes
        conns = region.domain.mesh.conns # connection of nodes in elements

        # define strings for evaluation
        funstr = '' # string to evaluate integral from variable values
        if self.desc.fun_args == 'einsum':
            funstr = "np." + self.desc.fun_args + "("
            from numpy import einsum
            funstr += "%s" % (self.desc.fun_summation)
        else:
            funstr = self.desc.fun_args + "("

        funvarstr = [] # string to evaluate variable at quadrature points
        indvar = 0 # index of variable argument
        indpar = 0 # index of parameter variable
        indlin = 0 # index of linear variable
        for ii, arg_name in enumerate(self.arg_names):
            if arg_name in vars_lin_names: # linear argument
                funvarstr.append('args[%d].field.eval_basis(maps)' \
                                 % (ii,))
                funstr += ", var_datas[%d][iqp, ilin[%d]]" % (ii, indlin)
                indvar += 1
                indlin += 1
            elif arg_name in self.names.material_name: # material argument
                mat_key = (self.region.name, self.integral_name)
                pom_str = "args[%d].datas[%s][maps.ig]['%s']" \
                    % (ii, str(mat_key), self.arg_methods[ii])
                pom_str += "[maps.iel].squeeze()"
                funvarstr.append(pom_str)
                funstr += ", var_datas[%d][iqp]" % (ii,)
            else: # parameter argument
                pom_str = 'np.tensordot(vec_el[%d], args[%d]' % (indpar, ii)
                pom_str += '.field.eval_basis(maps), (0, 1))'
                funvarstr.append(pom_str)
                funstr += ", var_datas[%d][iqp]" % (ii,)
                indvar += 1
                indpar += 1
        funstr += ")"
        fun = eval("lambda var_datas, iqp, ilin:" + funstr)

        funvar = []
        for ii in range(n_args):
            funvar.append(eval('lambda vec_el, args, maps:' + funvarstr[ii]))

        # loop over element groups
        for ig in self.iter_groups():
            gel_name = region.domain.groups[ig].gel.name
            qp_coors, qp_weights = self.integral.get_qp(gel_name)
            maps = Mapping(coor=np.take(coors, conns[ig][0], axis=0),
                           qp_coors=qp_coors, ig=ig, values=None)

            # creating local matrix
            val_loc_shape = []
            for arg in vars_lin:
                val_loc_shape.append(arg.field.aps[ig].econn.shape[1])
            if len(val_loc_shape) == 0:
                val_loc_shape = 1

            # material data
            mat_datas = []
            for mat in arg_mats:
                mat_datas.append(mat.datas)

            # set basis functions at reference element
            for ii, var in enumerate(arg_vars):
                var.field.set_basis(maps, self.arg_vars_method[ii])

            # store connectivity linear variables
            econn_lin = []
            for var in vars_lin:
                econn_lin.append(var.field.get_full_econn(ig))

            econn_par = []
            for var in vars_par:
                econn_par.append(var.field.get_full_econn(ig))

            # define iterable over basis functions of linear arguments
            if n_args_lin > 0:
                ilin_iter = np.arange(vars_lin[0].field.n_basis)[np.newaxis]
                for ii in range(1, n_args_lin):
                    n_bas = vars_lin[ii].field.n_basis
                    ilin_ii = np.tile(np.arange(n_bas), ilin_iter.shape[1])
                    ilin_iter = np.vstack([np.repeat(ilin_iter, n_bas),
                                          ilin_ii])
                ilin_iter = ilin_iter.T
            else:
                ilin_iter = [0]

            # loop over elements
            n_el = conns[ig].shape[0]
            for iel in np.arange(n_el):
                maps.update(np.take(coors, conns[ig][iel], axis=0), iel=iel)

                # get local vector at element from parameter variables
                vec_el = []
                for ii, arg in enumerate(vars_par):
                    vec_el.append(np.take(arg.data[0][arg.indx],
                                          econn_par[ii][iel]))

                # eval basis function in actual element
                var_datas = []
                for ii in range(n_args):
                    var_datas.append(funvar[ii](vec_el, args, maps))

                # local connectivity of linear variables
                econn_el = []
                for ii in np.arange(n_args_lin):
                    econn_el.append(econn_lin[ii][iel].tolist())

                # evaluate and assemble local matrix
                # loop over basis functions of linear variables
                for ilin in ilin_iter:
                    val_loc = 0.
                    for iqp in np.arange(maps.n_qp):
                        coef = qp_weights[iqp]*maps.det_jac
                        val_loc += coef*fun(var_datas, iqp, ilin)
                    # assemble to global matrix
                    asm_fun(asm_obj_T, econn_el, val_loc)

        if active_dof:
            eqi = self.get_eqi(kwargs, vars_lin_names)
            if n_args_lin == 1:
                asm_obj_T = asm_obj_T[eqi[0]]
            elif n_args_lin == 2:
                asm_obj_T = asm_obj_T.tocsr()[eqi[0][:, np.newaxis], eqi[1]]

        if n_args_lin>0:
            slc_global, n_global = self.get_slc_global(kwargs, vars_lin_names,
                                                       active_dof)

        if n_args_lin == 0:
            asm_obj_G = asm_obj_T
        elif n_args_lin == 1:
            asm_obj_G = np.zeros(n_global)
            asm_obj_G[slc_global[0]] = asm_obj_T
        elif n_args_lin == 2:
            asm_obj_G = np.zeros((n_global, n_global))
            asm_obj_G[slc_global[0], slc_global[1]] = asm_obj_T.todense()
            asm_obj_G = csr_matrix(asm_obj_G)

        # assemble linear object (asm_obj_G) to predefined one (asm_obj)
        if asm_obj is not None:
            if n_args_lin == 2:
                asm_obj = asm_obj_G
            else:
                asm_obj.data = asm_obj_G.data
        else:
            asm_obj = asm_obj_G

        return asm_obj

    @staticmethod
    def get_eqi(kwargs, vars_lin_names):
        eqi = []
        for var_name in vars_lin_names:
            var = kwargs[var_name]
            if var.is_state():
                eqi.append(var.eq_map.eqi)
            elif var.is_virtual():
                svar = var.get_primary()
                eqi.append(svar.eq_map.eqi)
            else:
                msg = "The variable '%s' is not supported as linear variable \
                    in IntFE term." % (var_name,)
                ValueError(msg)
        return eqi

    @staticmethod
    def get_slc_global(kwargs, vars_lin_names, active_dof):
        slc_global = []
        variables = kwargs[vars_lin_names[0]]._variables

        if active_dof:
            di = variables.adi
        else:
            di = variables.di

        for var_name in vars_lin_names:
            var = kwargs[var_name]
            if var.is_virtual():
                var_name = var.primary_var_name
            slc_global.append(di.indx[var_name])

        n_global = di.ptr[-1]

        return slc_global, n_global

    @staticmethod
    def get_dofs(kwargs, vars_lin_names):
        n_dof = []
        n_adof = []
        for var_name in vars_lin_names:
            var = kwargs[var_name]
            if var.is_state():
                pass
            elif var.is_virtual():
                var = var.get_primary()
            else:
                msg = "The variable '%s' is not supported as linear variable \
                    in IntFE term." % (var_name,)
                ValueError(msg)
            n_dof.append(var.n_dof)
            n_adof.append(var.n_dof)
        return n_dof, n_adof

    def get_conn_info(self):
        from sfepy.terms.terms import ConnInfo

        vvar = self.get_virtual_variable()
        svars = self.get_state_variables()
        pvars = self.get_parameter_variables()

        all_vars = self.get_variables()

        dc_type = self.get_dof_conn_type()
        tgs = self.get_geometry_types()

        v_igs = v_tg = None
        if vvar is not None:
            field = vvar.get_field()
            if field is not None:
                v_igs = field.igs
                if vvar.name in tgs:
                    v_tg = tgs[vvar.name]

                else:
                    v_tg = None

        else:
            # No virtual variable -> all unknowns are in fact known parameters.
            pvars += svars
            svars = []

        region = self.get_region()

        vals = []
        aux_pvars = []
        for svar in svars:
            # Allow only true state variables.
            if not svar.is_state():
                aux_pvars.append(svar)
                continue

            field = svar.get_field()
            if field is not None:
                s_igs = field.igs
            else:
                s_igs = None

            is_trace = False

            if svar.name in tgs:
                ps_tg = tgs[svar.name]
            else:
                ps_tg = v_tg

            val = ConnInfo(virtual=vvar, virtual_igs=v_igs,
                           state=svar, state_igs=s_igs,
                           primary=svar, primary_igs=s_igs,
                           has_virtual=True,
                           has_state=True,
                           is_trace=is_trace,
                           dc_type=dc_type,
                           v_tg=v_tg,
                           ps_tg=ps_tg,
                           region=region,
                           all_vars=all_vars)
            vals.append(val)

        pvars += aux_pvars
        for pvar in pvars:
            field = pvar.get_field()
            if field is not None:
                p_igs = field.igs
            else:
                p_igs = None
            is_trace = self.arg_traces[pvar.name]

            if pvar.name in tgs:
                ps_tg = tgs[pvar.name]
            else:
                ps_tg = v_tg

            val = ConnInfo(virtual=vvar, virtual_igs=v_igs,
                           state=None, state_igs=[],
                           primary=pvar.get_primary(), primary_igs=p_igs,
                           has_virtual=vvar is not None,
                           has_state=False,
                           is_trace=is_trace,
                           dc_type=dc_type,
                           v_tg=v_tg,
                           ps_tg=ps_tg,
                           region=region,
                           all_vars=all_vars)
            vals.append(val)

        if vvar and (len(vals) == 0):
            # No state, parameter variables, just the virtual one.
            val = ConnInfo(virtual=vvar, virtual_igs=v_igs,
                           state=vvar.get_primary(), state_igs=v_igs,
                           primary=vvar.get_primary(), primary_igs=v_igs,
                           has_virtual=True,
                           has_state=False,
                           is_trace=False,
                           dc_type=dc_type,
                           v_tg=v_tg,
                           ps_tg=v_tg,
                           region=region,
                           all_vars=all_vars)
            vals.append(val)

        return vals

class Mapping():
    """
    Class storing values useful for mapping between
    reference and true element.

    Values
    ------
    self.coor : numpy.array
        Coordinates of true element.
    self.qp_coors : numpy.array
        Coordinates of quadrature points.
    self.ig : int
        Actual index of group elements.
    self.jac : numpy.array
        Jacobian matrix of a mapping between reference and true element.
    self.det_jac : float
        Determinant of self.jac
    self.volume : float
        Volume of a true element.
    """
    def __init__(self, coor=None, qp_coors=None, ig=None, values=None):
        self.coor = coor
        self.ig = ig
        if qp_coors is None:
            self.qp_coor = coor
        else:
            self.qp_coor = qp_coors
        self.n_qp = self.qp_coor.shape[0]
        self.dim = coor.shape[1]
        self.n_coor = coor.shape[0]
        self.jac = np.zeros([self.dim, self.dim])
        self.det_jac = 0
        self.volume = 0

    def update(self, coor, iel=None):
        """
        Updates the values of reference mapping for particular element
        determined by coordinates.

        Input
        -----
        coor : numpy.array
            Coordinates of true element.
        """
        self.iel = iel
        self.get_Jacobian(coor)
        self.inv_jac = np.linalg.inv(self.jac)
        self.det_jac = np.linalg.det(self.jac)
        self.volume = 2*self.det_jac

    def get_Jacobian(self, coor):
        for ii in np.arange(self.dim):
            self.jac[:, ii] = coor[ii+1] - coor[0]

    def __repr__(self):
        ss = "Class : %s\n" % (self.__class__.__name__)
        ss += '    dim = %f \n' %(self.dim)
        ss += '    n_qp = %f \n' %(self.n_qp)
        ss += '    volume = %f \n' %(self.volume)
        ss += '    determinant of jacobian = %f \n' %(self.det_jac)
        return ss

def proceed_methods(bfref, methods):
    """
    It proceeds methods to array of basis functions evaluated at quadrature
    points.

    Parameters
    ----------
    bfref : numpy.array of shape = (n_qp, n_basis) + basis_shape
        An array that stores basis functions evaluated at quadrature points.
        Here n_qp denotes number of quadrature points, n_basis number of basis
        functions, and basis_shape is a shape of basis function, i.e.
        (1,) for scalar-valued,
        (dim,) for vector-valued,
        (dim, dim) for matrix-valued, etc.
    """

    for method in methods:

        if method in ['sym', 'Sym']:
            if bfref.ndim == 4 and bfref.shape[2] == bfref.shape[2]:
                bfref = 0.5*(bfref + bfref.swapaxes(-2, -1))
            else:
                raise ValueError("Improper shape of bfref")

        elif method in ['Man', 'Mandel', 'mandel']:
            from sfepy.mechanics.matcoefs import ElasticTensor
            bfref0 = np.copy(bfref)
            try:
                res_shape = ElasticTensor.create_mandel(bfref0[0, 0]).shape
                bfref = np.zeros(bfref0.shape[0:2] + res_shape)
            except:
                raise ValueError("Improper shape to create Mandel's notation.")

            for ii in np.arange(bfref.shape[0]):
                for jj in np.arange(bfref.shape[1]):
                    bfref[ii, jj] = ElasticTensor.create_mandel(bfref0[ii, jj])
        else:
            msg = "Improper method (%s) to proceed basis functions." \
                % str(method)
            raise ValueError(msg)

    return bfref

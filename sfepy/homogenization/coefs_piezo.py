from sfepy.base.base import *
from sfepy.homogenization.coefs_base import CoefSymSym, CoefDimSym, CoefDimDim, \
     CorrDimDim, CorrDim, ShapeDimDim, ShapeDim

class CorrectorsPiezoRS( CorrDimDim ):

    def get_variables( self, ir, ic, data ):
        """data: pis"""
        pis = data[self.requires[0]]
        yield (self.variables[-1], pis.states[ir,ic])

class CorrectorsPiezoK( CorrDim ):

    def get_variables( self, ir, data ):
        """data: pis_scalar"""
        pis_scalar = data[self.requires[0]]
        yield (self.variables[-1], pis_scalar.states[ir])

class ElasticPiezoCoef( CoefSymSym ):
    """Homogenized elastic tensor $E_{ijkl}$ for piezo-materials."""

    mode2var = {'row' : 0, 'col' : 1}

    def get_variables( self, problem, ir, ic, data, mode ):

        pis, corrs = [data[ii] for ii in self.requires]

        var_name = self.variables[self.mode2var[mode]]
        c_name = problem.variables[var_name].primary_var_name

        indx = corrs.di.indx[c_name]
        omega = corrs.states[ir,ic][indx]
        pi = pis.states[ir,ic] + omega
        yield (var_name, pi)

        var_name = self.variables[self.mode2var[mode] + 2]
        c_name = problem.variables[var_name].primary_var_name

        indx = corrs.di.indx[c_name]
        pi = corrs.states[ir,ic][indx]
        yield (var_name, pi)

class DielectricCoef( CoefDimDim ):

    def get_variables( self, problem, ir, ic, data, mode ):

        pis, corrs = [data[ii] for ii in self.requires]

        if mode == 'row':
            m2v = 0
            ii = ir
        else:
            m2v = 1
            ii = ic

        var_name = self.variables[m2v]
        c_name = problem.variables[var_name].primary_var_name

        indx = corrs.di.indx[c_name]
        pi = corrs.states[ii][indx]
        yield (var_name, pi)

        var_name = self.variables[m2v + 2]
        c_name = problem.variables[var_name].primary_var_name

        indx = corrs.di.indx[c_name]
        omega = corrs.states[ii][indx]
        pi = pis.states[ii] + omega
        yield (var_name, pi)

class PiezoCouplingCoef( CoefDimSym ):

    def get_variables( self, problem, ir, ic, data, mode ):

        pis_scalar, pis, corrs = [data[ii] for ii in self.requires]

        if mode == 'row':
            var_name = self.variables[2]
            yield (var_name, pis_scalar[ir])

        else:
            var_name = self.variables[0]
            c_name = problem.variables[var_name].primary_var_name

            indx = corrs.di.indx[c_name]
            omega = corrs.states[ir,ic][indx]
            pi = pis.states[ir,ic] + omega
            yield (var_name, pi)

            var_name = self.variables[1]
            c_name = problem.variables[var_name].primary_var_name

            indx = corrs.di.indx[c_name]
            pi = corrs.states[ir,ic][indx]
            yield (var_name, pi)

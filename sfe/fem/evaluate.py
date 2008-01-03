from sfe.base.base import *
import extmods.fem as fem
from sfe.terms import Term
from region import Region
from equations import Equation, buildArgs
from integrals import Integrals, quadratures

##
# 02.10.2007, c
class Evaluator( Struct ):
    pass

##
# 02.10.2007, c
class BasicEvaluator( Evaluator ):
    ##
    # 02.10.2007, c
    def __init__( self, problem, mtx = None, **kwargs ):
        Evaluator.__init__( self, problem = problem, data = kwargs )
        if mtx is None:
            self.mtx = problem.mtxA
        else:
            self.mtx = mtx

    ##
    # 02.12.2005, c
    # 25.07.2006
    # 02.10.2007
    def evalResidual( self, vec ):
        status = 0
        try:
            pb = self.problem
            vecR = evalResiduals( vec, pb.equations, pb.conf.fe.chunkSize,
                                  **self.data )
        except StopIteration, exc:
            vecR = None
            status = exc.args[0]
            print 'error %d in term "%s" of equation "%s"!'\
                  % (status, exc.args[1].name, exc.args[2].desc)

        return vecR, status
            
    ##
    # 02.12.2005, c
    # 25.07.2006
    # 02.10.2007
    def evalTangentMatrix( self, vec, mtx = None ):
        status = 0
        try:
            pb = self.problem
            if mtx is None:
                mtx = self.mtx
            mtx.data[:] = 0.0
            mtx = evalTangentMatrices( vec, mtx,
                                       pb.equations, pb.conf.fe.chunkSize,
                                       **self.data )
        except StopIteration, exc:
            status = exc.args[0]
            print ('error %d in term "%s" of derivative of equation "%s"'
                   + ' with respect to variable "%s"!')\
                  % (status, exc.args[1].name, exc.args[2].desc, exc.args[3] )

        return mtx, status

    ##
    # 02.12.2005, c
    # 09.12.2005
    # 25.07.2006
    # 02.10.2007
    def updateVec( self, vec, delta ):
        self.problem.updateVec( vec, delta )

##
# 04.10.2007, c
class LCBCEvaluator( BasicEvaluator ):

    ##
    # 04.10.2007, c
    def __init__( self, problem, mtx = None, **kwargs ):
        BasicEvaluator.__init__( self, problem, mtx, **kwargs )
        self.opLCBC = problem.variables.getLCBCOperator()

    ##
    # 04.10.2007, c
    def evalResidual( self, vec ):
        vecR, status = BasicEvaluator.evalResidual( self, vec )
        vecRR = self.opLCBC.T * vecR
        return vecRR, status
            
    ##
    # 04.10.2007, c
    def evalTangentMatrix( self, vec, mtx = None ):
        mtx, status = BasicEvaluator.evalTangentMatrix( self, vec, mtx )
        mtxR = self.opLCBC.T * mtx * self.opLCBC
        mtxR = mtxR.tocsr()
        mtxR.sort_indices()
##         import pylab
##         from sfe.base.plotutils import spy
##         spy( mtxR )
##         pylab.show()
        print mtxR.__repr__()
        return mtxR, status

    ##
    # 04.10.2007, c
    def updateVec( self, vec, deltaR ):
        delta = self.opLCBC * deltaR
        BasicEvaluator.updateVec( self, vec, delta )

##
# 03.09.2007, c
def assembleVector( vec, equation, variables, materials,
                    chunkSize = 1000, **kwargs ):
    getADofConn = variables.getADofConn

    for term in equation.terms:
#        print '>>>>>>', term.name, term.sign
        varNames = term.getVariableNames()
        args = buildArgs( term, variables, materials, **kwargs )
        vn = term.getVirtualName()
        dcType = term.getDofConnType()
##         print args
##         print vn
        for regionName, ig in term.iterGroups():
#            print term.region.name, regionName, ig
            term.setCurrentGroup( regionName, ig )
            dc = getADofConn( vn, True, dcType, regionName, ig )
#            print vn, dc.shape
#            pause()
            for vecInEls, iels, status in term( chunkSize = chunkSize,
                                                **args ):
                if status != 0:
                    raise StopIteration( status, term, equation )
                fem.assembleVector( vec, vecInEls, iels, term.sign, dc )

##
# 03.09.2007, c
def assembleMatrix( mtx, equation, variables, materials,
                    chunkSize = 1000, groupCanFail = True, **kwargs ):
    if not sp.isspmatrix_csr( mtx ):
        raise TypeError, 'must be CSR matrix!'
    tmd = (mtx.data, mtx.indptr, mtx.indices)

    getADofConn = variables.getADofConn

    for term in equation.terms:
#        print '>>>>>>', term.name, term.sign
        varNames = term.getVariableNames()
        args = buildArgs( term, variables, materials, **kwargs )
        vn = term.getVirtualName()
        sns = term.getStateNames()
        dcType = term.getDofConnType()

        for regionName, ig in term.iterGroups():
            term.setCurrentGroup( regionName, ig )
            rdc = getADofConn( vn, True, dcType, regionName, ig )
#            print vn, rdc.shape
            for sn in sns:
                cdc = getADofConn( sn, False, dcType, regionName, ig )
#                print sn, cdc.shape
#                pause()
                for mtxInEls, iels, status in term( diffVar = sn,
                                                    chunkSize = chunkSize,
                                                    **args ):
                    if status != 0:
                        raise StopIteration( status, term, equation, varNameCol )
                    fem.assembleMatrix( tmd[0], tmd[1], tmd[2], mtxInEls,
                                        iels, term.sign, rdc, cdc )

##
# 01.10.2007, c
def evalTermOP( state, termDesc, problem, **kwargs ):
    """Convenience wrapper of evalTerm() in a context of ProblemDefinition
    instance."""
    return evalTerm( state, termDesc, problem.conf,
                     problem.domain, problem.variables, problem.materials,
                     chunkSize = problem.domain.shape.nEl, **kwargs )

##
# created:       03.01.2006
# last revision: 14.12.2007
def evalTerm( state, termDesc, conf, domain, variables, materials,
              funmod = None, chunkSize = 1000, termPrefixes = None,
              caches = None, retCaches = False,
              override = True, newGeometries = True,
              dwMode = 'vector', tangentMatrix = None,
              **kwargs ):
    """Evaluate a term. May not succeed!"""
    if termPrefixes is None:
        termPrefixes = {}
    if caches is None:
        caches = {}

    equation = Equation.fromDesc( 'tmp', termDesc, termPrefixes )
    equation.parseTerms( domain.regions, caches )
    equation.setupTermArgs( variables, materials, kwargs )
    for cache in caches.itervalues():
        cache.setMode( override = override )

    if newGeometries:
        iNames = equation.getTermIntegralNames()
        integrals = Integrals.fromConf( conf.integrals, iNames )
        integrals.setQuadratures( quadratures )
        
        geometries = {}
        equation.describeGeometry( geometries, variables, integrals )

    variables.dataFromState( state )
    # itype according to the first term in termDesc!
    term = equation.terms[0]
    if term.itype == 'dw':

        variables.setupDofConns()
        if dwMode == 'vector':
            residual = variables.createStrippedStateVector()
            assembleVector( residual, equation, variables, materials,
                            chunkSize, groupCanFail = False, **kwargs )
            if variables.hasLCBC:
                opLCBC = variables.opLCBC
                residual = opLCBC.T * residual
            retVal = residual

        elif dwMode == 'matrix':
            if tangentMatrix is None:
                tangentMatrix = variables.createMatrixGraph()

            tangentMatrix.data[:] = 0.0
            assembleMatrix( tangentMatrix, equation, variables, materials,
                            chunkSize, groupCanFail = False, **kwargs )
            if variables.hasLCBC:
                opLCBC = variables.opLCBC
                tangentMatrix = opLCBC.T * tangentMatrix * opLCBC
                tangentMatrix = tangentMatrix.tocsr()
                tangentMatrix.sort_indices()
            retVal = tangentMatrix

        else:
            print dwMode
            raise ValueError

    elif term.itype == 'd':
        val = 0.0

        for term in equation.terms:
            args = buildArgs( term, variables, materials, **kwargs )
            for regionName, ig in term.iterGroups():
                for aux, iels, status in term( chunkSize = chunkSize,
                                               **args ):
                    val += term.sign * aux
            retVal = val

    elif (term.itype == 'de') or (term.itype == 'dq'):
        val = None

        for term in equation.terms:
            args = buildArgs( term, variables, materials, **kwargs )
            for regionName, ig in term.iterGroups():
                term.setCurrentGroup( regionName, ig )
                for aux, iels, status in term( chunkSize = chunkSize,
                                               **args ):
                    if val is None:
                        val = term.sign * aux
                    else:
                        val = nm.concatenate( (val, term.sign * aux), axis = 0 )
        retVal = val

    else:
        raise NotImplementedError, 'unknown term int. type: %s' % term.itype

    if retCaches:
        return retVal, caches
    else:
        return retVal

##
# 21.11.2005, c
# 27.11.2005
# 02.12.2005
# 09.12.2005
# 14.12.2005
# 21.03.2006
# 24.07.2006
# 25.07.2006
# 22.08.2006
# 11.10.2006
# 16.02.2007
# 03.09.2007
def evalResiduals( state, equations, chunkSize = 1000,
                   **kwargs ):

    variables = equations.variables
    materials = equations.materials

    variables.dataFromState( state )

    residual = variables.createStrippedStateVector()
    equations.invalidateTermCaches()

    for equation in equations:
        assembleVector( residual, equation, variables, materials,
                        chunkSize = chunkSize, **kwargs )

    return residual

##
# 22.11.2005, c
# 26.11.2005
# 02.12.2005
# 09.12.2005
# 14.12.2005
# 21.03.2006
# 24.07.2006
# 25.07.2006
# 22.08.2006
# 11.10.2006
# 16.02.2007
# 03.09.2007
def evalTangentMatrices( state, tangentMatrix, equations, chunkSize = 1000,
                         **kwargs ):

    variables = equations.variables
    materials = equations.materials

    variables.dataFromState( state )

    for equation in equations:
        assembleMatrix( tangentMatrix, equation, variables, materials,
                        chunkSize = chunkSize, **kwargs )

    return tangentMatrix

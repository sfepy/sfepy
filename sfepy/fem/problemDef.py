import os.path as op

from sfepy.base.base import *

from mesh import Mesh
from domain import Domain
from fields import Fields
from variables import Variables
from materials import Materials
from equations import Equations
from integrals import Integrals
import fea as fea
from sfepy.solvers.ts import TimeStepper
from sfepy.fem.evaluate import BasicEvaluator, LCBCEvaluator
from sfepy.solvers import Solver
from init_sfepy import install_dir

##
# 29.01.2006, c
class ProblemDefinition( Struct ):
    """
    Problem definition, the top-level class holding all data necessary to solve
    a problem.

    Contains: mesh, domain, materials, fields, variables, equations, solvers
    """
    ##
    # c: 29.01.2006, r: 09.07.2008
    def fromConf( conf,
                  initFields = True, initVariables = True, initEquations = True,
                  initSolvers = True ):

        mesh = Mesh.fromFile( conf.fileName_mesh )

        
        eldesc_dir = op.join( install_dir, 'eldesc' )
        domain = Domain.fromMesh( mesh, eldesc_dir )
        domain.setupGroups()
        domain.fixElementOrientation()
        domain.setupNeighbourLists()
#        print domain
##         for gr in domain.groups.itervalues():
##             print gr.shape

        domain.createRegions( conf.regions, conf.funmod )
##         print regions
##         pause()

        materials = Materials.fromConf( conf.materials )
        materials.setupRegions( domain.regions )
    ##     print materials
    ##     pause()

    ##     print fields
    ##     pause()

        obj = ProblemDefinition( conf = conf,
                                 domain = domain,
                                 materials = materials,
                                 eldesc_dir = eldesc_dir )

        if initFields:
            obj.setFields( conf.fields )

            if initVariables:
                obj.setVariables( conf.variables )

                if initEquations:
                    obj.setEquations( conf.equations )

        if initSolvers:
            obj.setSolvers( conf.solvers, conf.options )

        return obj
    fromConf = staticmethod( fromConf )

    ##
    # 18.04.2006, c
    def copy( self, **kwargs ):
        if 'share' in kwargs:
            share = kwargs['share']
            
        obj = ProblemDefinition()
        for key, val in self.__dict__.iteritems():
##             print key
            if key in share:
                obj.__dict__[key] = val
            else:
                obj.__dict__[key] = copy( val )
        return obj

    ##
    # c: 23.04.2007, r: 09.07.2008
    def setFields( self, conf_fields = None ):
        conf_fields = getDefault( conf_fields, self.conf.fields )
        fields = Fields.fromConf( conf_fields )
        fields.readInterpolants( self.eldesc_dir )
        fields.setupApproximations( self.domain )
##         print fields
##         print fields[0].aps
##         pause()
        fields.setupGlobalBase()
        fields.setupCoors()
        self.fields = fields

##         self.saveFieldMeshes( '.' )
##         pause()

    ##
    # c: 26.07.2006, r: 14.04.2008
    def setVariables( self, conf_variables = None ):
        conf_variables = getDefault( conf_variables, self.conf.variables )
        variables = Variables.fromConf( conf_variables, self.fields )
        variables.setupDofInfo() # Call after fields.setupGlobalBase().
        self.variables = variables
        self.mtxA = None
        self.solvers = None
##         print variables.di
##         pause()
        
    
    ##
    # c: 18.04.2006, r: 13.06.2008
    def setEquations( self, conf_equations, user = None, cacheOverride = None,
                      keepSolvers = False ):
        equations = Equations.fromConf( conf_equations )
        equations.parseTerms( self.domain.regions )
        equations.setupTermArgs( self.variables, self.materials, user )

        iNames = equations.getTermIntegralNames()
        self.integrals = Integrals.fromConf( self.conf.integrals, iNames )
        self.integrals.setQuadratures( fea.collectQuadratures() )

        self.geometries = {}
        equations.describeGeometry( self.geometries, self.variables,
                                    self.integrals )

##         print self.geometries
##         pause()
        # Call after describeGeometry(), as it sets ap.surfaceData.
        self.variables.setupDofConns()

        if cacheOverride is None:
            if hasattr( self.conf.fe, 'cacheOverride' ):
                cacheOverride = self.conf.fe.cacheOverride
            else:
                cacheOverride = True
        equations.setCacheMode( cacheOverride )

        self.equations = equations

        if not keepSolvers:
            self.solvers = None

    ##
    # c: 16.10.2007, r: 20.02.2008
    def setSolvers( self, conf_solvers = None, options = None ):
        """If solvers are not set in options, use first suitable in
        conf_solvers."""
        conf_solvers = getDefault( conf_solvers, self.conf.solvers )
        self.solverConfs = {}
        for key, val in conf_solvers.iteritems():
            self.solverConfs[val.name] = val
        
        def _findSuitable( prefix ):
            for key, val in self.solverConfs.iteritems():
                if val.kind.find( prefix ) == 0:
                    return val
            return None

        def _getSolverConf( kind ):
            try:
                key = options[kind]
                conf = self.solverConfs[key]
            except:
                conf = _findSuitable( kind + '.' )
            return conf
        
        self.tsConf = _getSolverConf( 'ts' )
        self.nlsConf = _getSolverConf( 'nls' )
        self.lsConf = _getSolverConf( 'ls' )
        info =  'using solvers:'
        if self.tsConf:
            info += '\n                ts: %s' % self.tsConf.name
        if self.nlsConf:
            info += '\n               nls: %s' % self.nlsConf.name
        if self.lsConf:
            info += '\n                ls: %s' % self.lsConf.name
        if info != 'using solvers:':
            output( info )

    ##
    # Utility functions below.
    ##

    ##
    # 17.10.2007, c
    def getSolverConf( self, name ):
        return self.solverConfs[name]
    
    ##
    # 29.01.2006, c
    # 25.07.2006
    def createStateVector( self ):
        return self.variables.createStateVector()

    ##
    # c: 08.08.2006, r: 14.04.2008
    def updateBC( self, ts, conf_ebc, conf_epbc, conf_lcbc, funmod,
                  createMatrix = False ):
        """Assumes same EBC/EPBC/LCBC nodes for all time steps. Otherwise set
        createMatrix to True"""
        self.variables.equationMapping( conf_ebc, conf_epbc,
                                        self.domain.regions, ts, funmod )
        self.variables.setupLCBCOperators( conf_lcbc, self.domain.regions )
                
        self.variables.setupADofConns()
        if (self.mtxA is None) or createMatrix:
            self.mtxA = self.variables.createMatrixGraph()
##             import sfepy.base.plotutils as plu
##             plu.spy( self.mtxA )
##             plu.pylab.show()

    ##
    # c: 13.06.2008, r: 13.06.2008
    def getDefaultTS( self, t0 = None, t1 = None, dt = None, nStep = None,
                      step = None ):
        t0 = getDefault( t0, 0.0 )
        t1 = getDefault( t1, 1.0 )
        dt = getDefault( dt, 1.0 )
        nStep = getDefault( nStep, 1 )
        ts = TimeStepper( t0, t1, dt, nStep )
        ts.setStep( step )
        return ts

    ##
    # c: 22.02.2008, r: 13.06.2008
    def updateMaterials( self, ts = None, funmod = None, extraMatArgs = None ):
        if ts is None:
            ts = self.getDefaultTS( step = 0 )
        funmod = getDefault( funmod, self.conf.funmod )
        self.materials.timeUpdate( ts, funmod, self.domain, extraMatArgs )

    ##
    # c: 12.03.2007, r: 13.06.2008
    def timeUpdate( self, ts = None,
                    conf_ebc = None, conf_epbc = None, conf_lcbc = None,
                    funmod = None, createMatrix = False, extraMatArgs = None ):
        if ts is None:
            ts = self.getDefaultTS( step = 0 )
            
        conf_ebc = getDefault( conf_ebc, self.conf.ebcs )
        conf_epbc = getDefault( conf_epbc, self.conf.epbcs )
        conf_lcbc = getDefault( conf_lcbc, self.conf.lcbcs )
        funmod = getDefault( funmod, self.conf.funmod )
        self.updateBC( ts, conf_ebc, conf_epbc, conf_lcbc, funmod, createMatrix )
        self.updateMaterials( ts, funmod, extraMatArgs )

    ##
    # 29.01.2006, c
    # 25.07.2006
    # 19.09.2006
    def applyEBC( self, vec, forceValues = None ):
        self.variables.applyEBC( vec, forceValues )

    ##
    # 25.07.2006, c
    def updateVec( self, vec, delta ):
        self.variables.updateVec( vec, delta )
        
    ##
    # c: 18.04.2006, r: 07.05.2008
    def stateToOutput( self, vec, fillValue = None, varInfo = None,
                       extend = True ):
        """
        Transforms state vector 'vec' to an output dictionary, that can be
        passed as 'out' kwarg to Mesh.write(). 'vec' must have full size,
        i.e. all fixed or periodic values must be included.

        Example:
        >>> out = problem.stateToOutput( state )
        >>> problem.saveState( 'file.vtk', out = out )

        Then the  dictionary entries a formed by components of the state vector
        corresponding to the unknown variables, each transformed to shape
        (n_mesh_nod, n_dof per node) - all values in extra nodes are removed.
        """
        return self.variables.stateToOutput( vec, fillValue, varInfo, extend )

    ##
    # 26.07.2006, c
    # 22.08.2006
    def getMeshCoors( self ):
        return self.domain.getMeshCoors()

    ##
    # created: 26.07.2006
    # last revision: 21.12.2007
    def setMeshCoors( self, coors, updateState = False ):
        fea.setMeshCoors( self.domain, self.fields, self.geometries,
                          coors, updateState )

    ##
    # c: 02.04.2008, r: 02.04.2008
    def initTime( self, ts ):
        self.equations.initTime( ts )

    ##
    # 08.06.2007, c
    def advance( self, ts ):
        self.equations.advance( ts )

    ##
    # c: 01.03.2007, r: 23.06.2008
    def saveState( self, fileName, state = None, out = None,
                   fillValue = None, postProcessHook = None,
                   filePerVar = False, **kwargs ):
        extend = not filePerVar
        if (out is None) and (state is not None):
            out = self.stateToOutput( state,
                                      fillValue = fillValue, extend = extend )
            if postProcessHook is not None:
                out = postProcessHook( out, self, state, extend = extend )

        if filePerVar:
            import os.path as op

            meshes = {}
            for var in self.variables.iterState():
                rname = var.field.region.name
                if meshes.has_key( rname ):
                    mesh = meshes[rname]
                else:
                    mesh = Mesh.fromRegion( var.field.region, self.domain.mesh,
                                            localize = True )
                    meshes[rname] = mesh
                vout = {}
                for key, val in out.iteritems():
                    if val.varName == var.name:
                        vout[key] = val
                base, suffix = op.splitext( fileName )
                mesh.write( base + '_' + var.name + suffix,
                            io = 'auto', out = vout, **kwargs )
        else:
            self.domain.mesh.write( fileName, io = 'auto', out = out, **kwargs )

    ##
    # c: 19.09.2006, r: 27.02.2008
    def saveEBC( self, fileName, force = True, default = 0.0 ):
        output( 'saving ebc...' )
        state = self.createStateVector()
        state.fill( default )
        if force:
            vals = dictFromKeysInit( [self.variables.names[ii]
                                      for ii in self.variables.state] )
            for ii, key in enumerate( vals.iterkeys() ):
                vals[key] = ii + 1
            self.applyEBC( state, forceValues = vals )
        else:
            self.applyEBC( state )
        self.saveState( fileName, state, fillValue = default )
        output( '...done' )

    ##
    # created:       30.03.2007
    # last revision: 27.02.2008
    def saveRegions( self, fileNameTrunk ):

        output( 'saving regions...' )
        for region in self.domain.regions:
            output( region.name )
            aux = Mesh.fromRegion( region, self.domain.mesh, self.domain.ed,
                                   self.domain.fa )
            aux.write( '%s_%s.mesh' % (fileNameTrunk, region.name), io = 'auto' )
        output( '...done' )

    ##
    # created:       02.01.2008
    # last revision: 27.02.2008
    def saveRegionFieldMeshes( self, fileNameTrunk ):

        output( 'saving regions of fields...' )
        for field in self.fields:
            fregion = self.domain.regions[field.regionName]
            output( 'field %s: saving regions...' % field.name )

            for region in self.domain.regions:
                if not fregion.contains( region ): continue
                output( region.name )
                aux = Mesh.fromRegionAndField( region, field )
                aux.write( '%s_%s_%s.mesh' % (fileNameTrunk,
                                              region.name, field.name),
                           io = 'auto' )
            output( '...done' )
        output( '...done' )

    ##
    # c: 03.07.2007, r: 27.02.2008
    def saveFieldMeshes( self, fileNameTrunk ):

        output( 'saving field meshes...' )
        for field in self.fields:
            output( field.name )
            field.writeMesh( fileNameTrunk + '_%s' )
        output( '...done' )

    ##
    # c: 17.01.2008, r: 11.04.2008
    def getEvaluator( self, fromNLS = False, **kwargs ):
        """
        Either create a new Evaluator instance (fromNLS == False),
        or return an existing instace, created in a preceding call to
        ProblemDefinition.initSolvers().
        """
        if fromNLS:
            solvers = self.getSolvers()
            try:
                ev = solvers.nls.evaluator
            except AttributeError:
                output( 'call ProblemDefinition.initSolvers() or'
                        ' set fromNLS to False!' )
                raise
        else:
            if self.variables.hasLCBC:
                ev = LCBCEvaluator( self, **kwargs )
            else:
                ev = BasicEvaluator( self, **kwargs )
        return ev

    ##
    # c: 04.04.2008, r: 04.04.2008
    def initSolvers( self, nlsStatus = None, lsConf = None, nlsConf = None,
                     mtx = None, **kwargs ):
        lsConf = getDefault( lsConf, self.lsConf, 'you must set linear solver!' )

        nlsConf = getDefault( nlsConf, self.nlsConf,
                              'you must set nonlinear solver!' )
        
        ev = self.getEvaluator( **kwargs )
        ls = Solver.anyFromConf( lsConf, mtx = mtx )
        nls = Solver.anyFromConf( nlsConf, evaluator = ev,
                                  linSolver = ls, status = nlsStatus )
        self.solvers = Struct( name = 'solvers', ls = ls, nls = nls )

    ##
    # c: 04.04.2008, r: 04.04.2008
    def getSolvers( self ):
        return getattr( self, 'solvers', None )

    ##
    # c: 04.04.2008, r: 04.04.2008
    def isLinear( self ):
        nlsConf = getDefault( None, self.nlsConf,
                              'you must set nonlinear solver!' )

        if nlsConf.problem == 'linear':
            return True
        else:
            return False

    ##
    # c: 13.06.2008, r: 13.06.2008
    def setLinear( self, isLinear ):
        nlsConf = getDefault( None, self.nlsConf,
                              'you must set nonlinear solver!' )
        if isLinear:
            nlsConf.problem = 'linear'
        else:
            nlsConf.problem = 'nonlinear'

    ##
    # c: 03.10.2007, r: 11.04.2008
    def solve( self, state0 = None, nlsStatus = None,
               lsConf = None, nlsConf = None, forceValues = None,
               **kwargs ):
        """Solve self.equations in current time step."""
        solvers = self.getSolvers()
        if solvers is None:
            self.initSolvers( nlsStatus, lsConf, nlsConf, **kwargs )
            solvers = self.getSolvers()
        else:
            if kwargs:
                ev = self.getEvaluator( fromNLS = True )
                ev.setTermArgs( **kwargs )
            
        if state0 is None:
            state = self.createStateVector()
        else:
            state = state0.copy()

        self.applyEBC( state, forceValues = forceValues )

        state = solvers.nls( state )

        return state

    ##
    # c: 06.02.2008, r: 04.04.2008
    def getTimeSolver( self, tsConf = None, **kwargs ):
        tsConf = getDefault( tsConf, self.tsConf,
                             'you must set time-stepping solver!' )
        
        return Solver.anyFromConf( tsConf, **kwargs )

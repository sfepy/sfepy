from sfe.base.base import *

from mesh import Mesh
from domain import Domain
from fields import Fields
from variables import Variables
from materials import Materials
from equations import Equations
from integrals import Integrals
import fea as fea
from sfe.base.ioutils import writeVTK
from sfe.solvers.ts import TimeStepper
from sfe.fem.evaluate import BasicEvaluator, LCBCEvaluator
from sfe.solvers import Solver

##
# 29.01.2006, c
class ProblemDefinition( Struct ):
    
    ##
    # 29.01.2006, c
    # 21.03.2006
    # 18.04.2006
    # 08.06.2006
    # 25.07.2006
    # 22.08.2006
    # 24.08.2006
    # 19.09.2006
    # 23.04.2007
    # 16.10.2007
    def fromConf( conf,
                  initFields = True, initVariables = True, initEquations = True,
                  initSolvers = True ):

        mesh = Mesh.fromFile( conf.fileName_mesh )

        domain = Domain.fromMesh( mesh, 'eldesc' )
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
                                 materials = materials )

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
            print key
            if key in share:
                obj.__dict__[key] = val
            else:
                obj.__dict__[key] = copy( val )
        return obj

    ##
    # 23.04.2007, c
    # 04.09.2007
    def setFields( self, conf_fields = None ):
        conf_fields = getDefault( conf_fields, self.conf.fields )
        fields = Fields.fromConf( conf_fields )
        fields.readInterpolants( 'eldesc' )
        fields.setupApproximations( self.domain )
##         print fields
##         print fields[0].aps
##         pause()
        fields.setupGlobalBase()
        fields.setupCoors()
        self.fields = fields

#        self.saveFieldMeshes( '.' )

    ##
    # 26.07.2006, c
    # 08.08.2006
    # 11.07.2007
    # 04.09.2007
    def setVariables( self, conf_variables = None ):
        conf_variables = getDefault( conf_variables, self.conf.variables )
        variables = Variables.fromConf( conf_variables, self.fields )
        variables.setupDofInfo() # Call after fields.setupGlobalBase().
        self.variables = variables
        self.mtxA = None
##         print variables.di
##         pause()
        
    
    ##
    # 18.04.2006, c
    # 25.07.2006
    # 24.08.2006
    # 10.10.2006
    # 22.02.2007
    # 23.04.2007
    # 11.07.2007
    def setEquations( self, conf_equations, user = None, cacheOverride = True ):
        equations = Equations.fromConf( conf_equations )
        equations.parseTerms( self.domain.regions )
        equations.setupTermArgs( self.variables, self.materials, user )

        iNames = equations.getTermIntegralNames()
        self.integrals = Integrals.fromConf( self.conf.integrals, iNames )
        self.integrals.setQuadratures( fea.collectQuadratures() )

        self.geometries = {}
        equations.describeGeometry( self.geometries, self.variables,
                                    self.integrals )

        # Call after describeGeometry(), as it sets ap.surfaceData.
        self.variables.setupDofConns()

        equations.setCacheMode( cacheOverride )

        self.equations = equations

    ##
    # 16.10.2007, c
    # 17.10.2007
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
        print 'using solvers:'
        if self.tsConf:
            print '           ts: %s' % self.tsConf.name
        if self.nlsConf:
            print '          nls: %s' % self.nlsConf.name
        if self.lsConf:
            print '           ls: %s' % self.lsConf.name

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
    # 08.08.2006, c
    # 16.10.2006
    # 12.03.2007
    # 03.10.2007
    def updateBC( self, ts, conf_ebc, conf_epbc, conf_lcbc, funmod ):
        """Assumes same EBC/EPBC/LCBC nodes for all time steps."""
        self.variables.equationMapping( conf_ebc, conf_epbc,
                                        self.domain.regions, ts, funmod )
        self.variables.setupLCBCOperators( conf_lcbc, self.domain.regions )
                
        self.variables.setupADofConns()
        if self.mtxA is None:
            self.mtxA = self.variables.createMatrixGraph()
##             import sfe.base.plotutils as plu
##             plu.spy( self.mtxA )
##             plu.pylab.show()


    ##
    # 12.03.2007, c
    # 11.04.2007
    # 03.10.2007
    # 29.10.2007
    def timeUpdate( self, ts = None,
                    conf_ebc = None, conf_epbc = None, conf_lcbc = None,
                    funmod = None, extraMatArgs = None ):
        if ts is None:
            ts = TimeStepper( 0.0, 1.0, 1.0, 1 )
            ts.setStep( 0 )
            
        conf_ebc = getDefault( conf_ebc, self.conf.ebc )
        conf_epbc = getDefault( conf_epbc, self.conf.epbc )
        conf_lcbc = getDefault( conf_lcbc, self.conf.lcbc )
        funmod = getDefault( funmod, self.conf.funmod )
        self.updateBC( ts, conf_ebc, conf_epbc, conf_lcbc, funmod )
        self.materials.timeUpdate( ts, funmod, self.domain, extraMatArgs )

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
    # 18.04.2006, c
    # 25.07.2006
    def stateToOutput( self, vec, fillValue = None ):
        return self.variables.stateToOutput( vec, fillValue )

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
    # 08.06.2007, c
    def advance( self, ts ):
        self.equations.advance( ts )

    ##
    # 01.03.2007, c
    # 04.06.2007
    def saveStateToVTK( self, fileName, state, fillValue = None ):
        out = self.stateToOutput( state, fillValue )
        fd = open( fileName, 'w' )
        writeVTK( fd, self.domain.mesh, out )
        fd.close()

    ##
    # 19.09.2006, c
    # 01.03.2007
    # 09.03.2007
    # 15.03.2007
    # 04.06.2007
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
        self.saveStateToVTK( fileName, state, fillValue = default )

    ##
    # created:       30.03.2007
    # last revision: 14.12.2007
    def saveRegions( self, fileNameTrunk ):

        for region in self.domain.regions:
            output( 'saving region %s...' % region.name )
            aux = Mesh.fromRegion( region, self.domain.mesh, self.domain.ed,
                                   self.domain.fa )
            aux.write( '%s_%s.mesh' % (fileNameTrunk, region.name) )

    ##
    # created:       02.01.2008
    # last revision: 02.01.2008
    def saveRegionFieldMeshes( self, fileNameTrunk ):

        for field in self.fields:
            fregion = self.domain.regions[field.regionName]
            output( 'saving regions of field %s:' % field.name )

            for region in self.domain.regions:
                if not fregion.contains( region ): continue
                output( '  saving region %s...' % region.name )
                aux = Mesh.fromRegionAndField( region, field )
                aux.write( '%s_%s_%s.mesh' % (fileNameTrunk,
                                              region.name, field.name) )

    ##
    # 03.07.2007, c
    def saveFieldMeshes( self, fileNameTrunk ):

        for field in self.fields:
            output( 'saving field %s...' % field.name )
            field.writeMesh( fileNameTrunk + '_%s' )

    ##
    # 03.10.2007, c
    # 10.10.2007
    # 17.10.2007
    def solve( self, state0 = None, data = None, nlsStatus = None,
               lsConf = None, nlsConf = None ):
        """Solve self.equations in current time step."""

        lsConf = getDefault( lsConf, self.lsConf )
        if lsConf is None:
            print 'you must set linear solver!'
            raise ValueError

        nlsConf = getDefault( nlsConf, self.nlsConf )
        if nlsConf is None:
            print 'you must set nonlinear solver!'
            raise ValueError
        
        if state0 is None:
            state = self.createStateVector()
        else:
            state = state0.copy()

        self.applyEBC( state )
        if self.variables.hasLCBC:
            ev = LCBCEvaluator( self, data = data )
        else:
            ev = BasicEvaluator( self, data = data )

        ls = Solver.anyFromConf( lsConf )
        nls = Solver.anyFromConf( nlsConf, evaluator = ev,
                                  linSolver = ls, status = nlsStatus )
        state = nls( state )

        return state

from sfe.base.base import *

import sfe.base.ioutils as io
from sfe.fem.problemDef import ProblemDefinition

required = ['fileName_mesh', 'field_[0-9]+|fields',
            'ebc_[0-9]+|ebcs', 'fe', 'equations',
            'region_[0-9]+|regions', 'variable_[0-9]+|variables',
            'material_[0-9]+|materials', 'integral_[0-9]+|integrals',
            'solver_[0-9]+|solvers']
other = ['functions', 'modules', 'epbc_[0-9]+|epbcs',
         'lcbc_[0-9]+|lcbcs', 'nbc_[0-9]+|nbcs', 'options']


##
# c: 03.07.2007, r: 02.01.2008
def saveOnly( conf, saveNames, problem = None ):
    """Save information available prior to setting equations and
    solving them."""
    if problem is None:
        problem = ProblemDefinition.fromConf( conf, initVariables = False )

    if saveNames.regions is not None:
        problem.saveRegions( saveNames.regions )

    if saveNames.fieldMeshes is not None:
        problem.saveFieldMeshes( saveNames.fieldMeshes )

    if saveNames.regionFieldMeshes is not None:
        problem.saveRegionFieldMeshes( saveNames.regionFieldMeshes )

    if saveNames.ebc is not None:
        if not hasattr( problem, 'variables' ):
            problem.setVariables( conf.variables )
        try:
            ts = TimeStepper.fromConf( conf.ts )
            ts.setStep( 0 )
        except:
            ts = None
        try:
            problem.variables.equationMapping( conf.ebc, conf.epbc,
                                               problem.regions, ts,
                                               conf.funmod )
        except Exception, e:
            output( 'cannot make equation mapping!' )
            output( 'reason: %s' % e )
        else:
            problem.saveEBC( saveNames.ebc )

##
# 20.03.2007, c
# 30.03.2007
# 28.05.2007
# 03.07.2007
# 18.07.2007
# 02.10.2007
# 03.10.2007
def solveStationary( conf, data = None, saveNames = None, nlsStatus = None ):

    if data is None:
        # Term-dependent data.
        data = {}
    problem = ProblemDefinition.fromConf( conf )

    problem.timeUpdate( None )

    if saveNames is not None:
        saveOnly( conf, saveNames, problem = problem )

    state = problem.solve( nlsStatus = nlsStatus )

    return problem, state, data


##
# c: 06.02.2008, r: 06.02.2008
def prepareSaveData( ts, conf, options ):
    if options.outputFileNameTrunk:
        ofnTrunk = options.outputFileNameTrunk
    else:
        ofnTrunk = io.getTrunk( conf.fileName_mesh ) + '_out'

    suffix = '.%%0%dd.vtk' % ts.nDigit

    try:
        saveSteps = conf.options.saveSteps
    except:
        saveSteps = -1
    if saveSteps == -1:
        saveSteps = ts.nStep

    isSave = nm.linspace( 0, ts.nStep - 1, saveSteps ).astype( nm.int32 )
    isSave = nm.unique1d( isSave )

    return ofnTrunk, suffix, isSave

##
# c: 06.02.2008, r: 06.02.2008
def timeStepFunction( ts, state0, problem, data ):
    problem.timeUpdate( ts )

    vh = problem.conf.options.variableHistory
    varNames0 = vh.values()
    varNames = vh.keys()

    if ts.step == 0:
        state = state0.copy()
        problem.applyEBC( state )
    else:
        problem.variables.nonStateDataFromState( varNames0, state0, varNames )
        state = problem.solve( state0 = state0, ts = ts, **data )

    problem.advance( ts )

    return state

##
# c: 12.01.2007, r: 08.02.2008
def solveDirect( conf, options ):
    """Generic (simple) problem solver."""
    if options.outputFileNameTrunk:
        ofnTrunk = options.outputFileNameTrunk
    else:
        ofnTrunk = io.getTrunk( conf.fileName_mesh )

    saveNames = Struct( ebc = None, regions = None, fieldMeshes = None,
                        regionFieldMeshes = None )
    if options.saveEBC:
        saveNames.ebc = ofnTrunk + '_ebc.vtk'
    if options.saveRegions:
        saveNames.regions = ofnTrunk + '_region'
    if options.saveFieldMeshes:
        saveNames.fieldMeshes = ofnTrunk + '_field'
    if options.saveRegionFieldMeshes:
        saveNames.regionFieldMeshes = ofnTrunk + '_region_field'

    if options.solveNot:
        saveOnly( conf, saveNames )
        return None, None, None

    if hasattr( conf.options, 'ts' ):
        ##
        # Time-dependent problem.
        data = {}
        pb = ProblemDefinition.fromConf( conf )
        timeSolver = pb.getTimeSolver( stepFun = timeStepFunction,
                                       stepArgs = (pb, data) )

        ofnTrunk, suffix, isSave = prepareSaveData( timeSolver.ts,
                                                    conf, options )
        
        state0 = pb.createStateVector()
        ii = 0
        for step, time, state in timeSolver( state0 ):
            
            if isSave[ii] == step:
                pb.saveState( ofnTrunk + suffix % step, state )

                ii += 1
    else:
        ##
        # Stationary problem.
        pb, state, data = solveStationary( conf, saveNames = saveNames )
        pb.saveState( ofnTrunk + '.vtk', state )

        if options.dump:
            import tables as pt
            import numarray as nar

            fd = pt.openFile( ofnTrunk + '.h5', mode = 'w', title = "Dump file" )
            for key, val in out.iteritems():
                fd.createArray( fd.root, key, nar.asarray( val.data ), 
                                '%s data' % val.mode )
            fd.close()

    return pb, state, data

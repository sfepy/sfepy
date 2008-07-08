from sfepy.base.base import *

import sfepy.base.ioutils as io
from sfepy.fem.problemDef import ProblemDefinition
from sfepy.base.conf import getStandardKeywords

##
# c: 03.07.2007, r: 27.02.2008
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
            problem.variables.equationMapping( conf.ebcs, conf.epbcs,
                                               problem.domain.regions, ts,
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
# c: 06.02.2008, r: 08.07.2008
def timeStepFunction( ts, state0, problem, data ):
    problem.timeUpdate( ts )

    vh = getDefaultAttr( problem.conf.options, 'variableHistory', {} )
    varNames0, varNames = vh.values(), vh.keys()

    setHistoryData = problem.variables.nonStateDataFromState

    if ts.step == 0:
        state = state0.copy()
        problem.initTime( ts )
        problem.applyEBC( state )

        if problem.equations.caches:
            # Initialize caches.
            ev = problem.getEvaluator( ts = ts, **data )
            setHistoryData( varNames0, state0, varNames )
            vecR, ret = ev.evalResidual( state )
            if ret == 0: # OK.
                err = nla.norm( vecR )
                output( 'initial residual: %e' % err )
            else:
                output( 'initial residual evaluation failed, giving up...' )
                raise ValueError
        if problem.isLinear():
            # Assemble linear system matrix for all
            # time steps.
            ev = problem.getEvaluator( ts = ts, mtx = problem.mtxA, **data )
            setHistoryData( varNames0, state0, varNames )
            mtxA, ret = ev.evalTangentMatrix( state )
            if ret != 0:
                output( 'matrix evaluation failed, giving up...' )
                raise ValueError
        else:
            mtxA = None

        # Initialize solvers (and possibly presolve the matrix).
        problem.initSolvers( ts = ts, mtx = mtxA, **data )

    else:
        setHistoryData( varNames0, state0, varNames )
        state = problem.solve( state0 = state0, ts = ts, **data )

    problem.advance( ts )

    return state

##
# c: 13.06.2008, r: 13.06.2008
def solveEvolutionaryOP( problem, options,
                         saveResults = True, returnHistory = False,
                         postProcessHook = None ):
    """TODO  returnHistory"""
    
    data = {}
    timeSolver = problem.getTimeSolver( stepFun = timeStepFunction,
                                        stepArgs = (problem, data) )

    ofnTrunk, suffix, isSave = prepareSaveData( timeSolver.ts,
                                                problem.conf, options )

    state0 = problem.createStateVector()
    ii = 0
    for step, time, state in timeSolver( state0 ):

        if saveResults and (isSave[ii] == step):
            problem.saveState( ofnTrunk + suffix % step, state,
                               postProcessHook = postProcessHook )

            ii += 1
    return state, data

##
# c: 13.06.2008, r: 13.06.2008
def solveStationaryOP( problem, options, saveResults = True, ts = None,
                       postProcessHook = None ):
    data = {}
    problem.timeUpdate( ts )
    state = problem.solve()

    if saveResults:
        if options.outputFileNameTrunk:
            ofnTrunk = options.outputFileNameTrunk
        else:
            ofnTrunk = io.getTrunk( problem.conf.fileName_mesh ) + '_out'
        problem.saveState( ofnTrunk + '.vtk', state,
                           postProcessHook = postProcessHook )

    return state, data
    
##
# c: 12.01.2007, r: 22.06.2008
def solveDirect( conf, options ):
    """Generic (simple) problem solver."""
    if options.outputFileNameTrunk:
        ofnTrunk = options.outputFileNameTrunk
    else:
        ofnTrunk = io.getTrunk( conf.fileName_mesh )

    opts = conf.options
    if hasattr( opts, 'postProcessHook' ) and opts.postProcessHook is not None:
        # User postprocessing.
        postProcessHook = getattr( conf.funmod, opts.postProcessHook )
    else:
        postProcessHook = None

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

    isExtraSave = False
    for name, val in saveNames.toDict().iteritems():
        if val is not None:
            isExtraSave = True
            break
    if isExtraSave:
        saveOnly( conf, saveNames )

    if options.solveNot:
        return None, None, None
            
    pb = ProblemDefinition.fromConf( conf )
    if hasattr( conf.options, 'ts' ):
        ##
        # Time-dependent problem.
        state, data = solveEvolutionaryOP( pb, options,
                                           postProcessHook = postProcessHook )
    else:
        ##
        # Stationary problem.
        state, data = solveStationaryOP( pb, options,
                                         postProcessHook = postProcessHook )

        if options.dump:
            import tables as pt
            import numarray as nar

            fd = pt.openFile( ofnTrunk + '.h5', mode = 'w', title = "Dump file" )
            for key, val in out.iteritems():
                fd.createArray( fd.root, key, nar.asarray( val.data ), 
                                '%s data' % val.mode )
            fd.close()

    return pb, state, data

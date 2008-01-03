from sfe.base.base import *
from sfe.fem.problemDef import ProblemDefinition
from sfe.solvers.ts import TimeStepper

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

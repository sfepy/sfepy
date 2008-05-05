# c: 05.05.2008, r: 05.05.2008
import sys
sys.path.append( '.' )

from sfe.fem.periodic import *

# c: 05.05.2008, r: 05.05.2008
def defineRegions( fileName ):
    """Define various subdomain for a given mesh file. This function is called
    below."""
    regions = {}
    is3D = False
    
    regions['Y'] = ('all', {})

    eog = 'elements of group %d'
    if fileName.find( 'osteonT1' ) >= 0:
        matIds = [11, 39, 6, 8, 27, 28, 9, 2, 4, 14, 12, 17, 45, 28, 15]
        regions['Ym'] = (' +e '.join( (eog % im) for im in  matIds ), {})
        wx = 0.865
        wy = 0.499

    regions['Yc'] = ('r.Y -e r.Ym', {})

    # Sides.
    regions['Left'] = ('nodes in (x < -%.3f)' % wx, {})
    regions['Right'] = ('nodes in (x > %.3f)' % wx, {})
    regions['Bottom'] = ('nodes in (y < -%.3f)' % wy, {})
    regions['Top'] = ('nodes in (y > %.3f)' % wy, {})
    regions['Corners'] = ("""nodes in
                            ((x < -%.3f) & (y < -%.3f))
                          | ((x >  %.3f) & (y < -%.3f))
                          | ((x >  %.3f) & (y >  %.3f))
                          | ((x < -%.3f) & (y >  %.3f))
                          """ % ((wx, wy) * 4), {})
    return is3D, regions, matIds

##
# c: 05.05.2008, r: 05.05.2008
def getPars( ts, coor, region, ig, matIds = [] ):
    """Define material parameters:
         $D_ijkl$ (elasticity),
       in a given region."""
    dim = coor.shape[1]
    sym = (dim + 1) * dim / 2

    m2i = region.domain.matIdsToIGs
    matrixIgs = [m2i[im] for im in matIds]

    out = {}

    # in 1e+10 [Pa]
    lam = 1.7
    mu = 0.3
    o = nm.array( [1.] * dim + [0.] * (sym - dim), dtype = nm.float64 )
    oot = nm.outer( o, o )
    out['D'] = lam * oot + mu * nm.diag( o + 1.0 )

    if ig not in matrixIgs: # channels
        out['D'] *= 1e-1

    return out
    
##
# Mesh file.
fileName_mesh = 'examples/osteonT1_11.mesh'

##
# Define regions (subdomains, boundaries) - $Y$, $Y_i$, ...
# depending on a mesh used.
is3D, regions, matIds = defineRegions( fileName_mesh )

if is3D:
    dim, geom = 3, '3_4'
else:
    dim, geom = 2, '2_3'

##
# Define fields: 'displacement' in $Y$,
# 'pressure_m' in $Y_m$.
field_1 = {
    'name' : 'displacement',
    'dim' : (dim,1),
    'domain' : 'Y',
    'bases' : {'Y' : '%s_P1' % geom}
}

field_2 = {
    'name' : 'pressure_m',
    'dim' : (1,1),
    'domain' : 'Ym',
    'bases' : {'Ym' : '%s_P1' % geom}
}

##
# Define corrector variables: unknown displaements: uc, test: vc
# displacement-like variables: Pi, Pi1, Pi2
variables = {
    'uc'       : ('unknown field',   'displacement', 0),
    'vc'       : ('test field',      'displacement', 'uc'),
    'Pi'       : ('parameter field', 'displacement', 'uc'),
    'Pi1'      : ('parameter field', 'displacement', 'uc'),
    'Pi2'      : ('parameter field', 'displacement', 'uc'),
}

##
# Periodic boundary conditions.
if dim == 3:
    epbc_10 = {
        'name' : 'periodic_x',
        'region' : ['Left', 'Right'],
        'dofs' : {'uc.all' : 'uc.all'},
        'match' : 'matchXPlane',
    }
    epbc_11 = {
        'name' : 'periodic_y',
        'region' : ['Near', 'Far'],
        'dofs' : {'uc.all' : 'uc.all'},
        'match' : 'matchYPlane',
    }
    epbc_12 = {
        'name' : 'periodic_z',
        'region' : ['Top', 'Bottom'],
        'dofs' : {'uc.all' : 'uc.all'},
        'match' : 'matchZPlane',
    }
else:
    epbc_10 = {
        'name' : 'periodic_x',
        'region' : ['Left', 'Right'],
        'dofs' : {'uc.all' : 'uc.all'},
        'match' : 'matchYLine',
    }
    epbc_11 = {
        'name' : 'periodic_y',
        'region' : ['Top', 'Bottom'],
        'dofs' : {'uc.all' : 'uc.all'},
        'match' : 'matchXLine',
    }
    
##
# Dirichlet boundary conditions.
ebcs = {
    'fixed_u' : ('Corners', {'uc.all' : 0.0}),
}

##
# Material defining constitutive parameters of the microproblem.
material_1 = {
    'name' : 'm',
    'mode' : 'function',
    'region' : 'Y',
    'function' : 'getPars',
    'extraArgs' : {'matIds' : matIds},
}

##
# Numerical quadratures for volume (i3 - order 3) integral terms.
integral_1 = {
    'name' : 'i3',
    'kind' : 'v',
    'quadrature' : 'gauss_o3_d%d' % dim,
}

##
# Steady state correctors $\bar{\omega}^{rs}$.
equations = {
    'eq_1' : 
    """dw_lin_elastic.i3.Y( m.D, vc, uc )
       = - dw_lin_elastic_r.i3.Y( m.D, vc, Pi )""",
}

##
# FE assembling options.
fe = {
    'chunkSize' : 100000,
    'cacheOverride' : True,
}

##
# Solvers.
solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.umfpack', # Direct solver.
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'iMax'      : 2,
    'epsA'      : 1e-8,
    'epsR'      : 1e-2,
    'macheps'   : 1e-16,
    'linRed'    : 1e-2, # Linear system error < (epsA * linRed).
    'lsRed'     : 0.1,
    'lsRedWarp' : 0.001,
    'lsOn'      : 0.99999,
    'lsMin'     : 1e-5,
    'check'     : 0,
    'delta'     : 1e-6,
    'isPlot'    : False,
    'problem'   : 'nonlinear', # 'nonlinear' or 'linear' (ignore iMax)
}

############################################
# Mini-application below, computing the homogenized elastic coefficients.

##
# c: 28.02.2007, r: 13.02.2008
def buildOpPi( varName, problem, ir, ic ):
    """\Pi^{rs}_i = \delta_{ri} y_s. """
    var = problem.variables[varName]
    coor = var.field.getCoor()[:,:-1]

    pi = nm.zeros_like( coor )
    pi[:,ir] = coor[:,ic]
    pi.shape = (pi.shape[0] * pi.shape[1],)

    return pi

##
# c: 05.05.2008, r: 05.05.2008
def createPis( problem, variables, varName ):
    problem.setVariables( variables )

    dim = problem.domain.mesh.dim
    pis = nm.zeros( (dim, dim), dtype = nm.object )
    for ir in range( dim ):
        for ic in range( dim ):
            pi = buildOpPi( varName, problem, ir, ic )
            pis[ir,ic] = pi
    return pis

##
# c: 05.05.2008, r: 05.05.2008
def  solveSteadyCorrectors_rs( problem, equations, variables, pis,
                               ofnTrunk, postProcessHook = None,
                               filePerVar = False ):
    """Compute the steady state correctors $\bar{\omega}^{rs}$"""
    from sfe.base.base import Struct
    
    dim = problem.domain.mesh.dim

    problem.setVariables( variables )
    problem.setEquations( equations )

    problem.timeUpdate()

    statesRS = nm.zeros( (dim, dim), dtype = nm.object )
    for ir in range( dim ):
        for ic in range( dim ):
            pi = pis[ir,ic]
            # Non-state variables must be assigned manually.
            problem.variables['Pi'].dataFromData( pi )

            state = problem.createStateVector()
            problem.applyEBC( state )
            state = problem.solve()
            assert problem.variables.hasEBC( state )
            statesRS[ir,ic] = state

            problem.saveState( ofnTrunk + '_steady_rs_%d%d.vtk' % (ir, ic),
                               state, postProcessHook = postProcessHook,
                               filePerVar = filePerVar )
    return Struct( name = 'Steady RS correctors',
                   statesRS = statesRS,
                   di = problem.variables.di )

##
# c: 05.03.2008, r: 05.03.2008
def iterSym( dim ):
    for ii in xrange( dim ):
        yield ii, ii
    for ir in xrange( 0, dim ):
        for ic in xrange( ir + 1, dim ):
            yield ir, ic

##
# c: 05.05.2008, r: 05.05.2008
def coefE( problem, corrsRS, pis ):
    """Homogenized elastic coefficient $E_{ijkl}$."""
    from sfe.fem.evaluate import evalTermOP

    coefTerm = 'd_lin_elastic.i3.Y( m.D, Pi1, Pi2 )'

    dim = problem.domain.mesh.dim
    sym = (dim + 1) * dim / 2
    coef = nm.zeros( (sym, sym), dtype = nm.float64 )

    indx = corrsRS.di.indx['uc']
    for ir, (irr, icr) in enumerate( iterSym( dim ) ):
        omega1 = corrsRS.statesRS[irr,icr][indx]
        pi1 = pis[irr,icr] + omega1
        # Non-state variables must be assigned manually.
        problem.variables['Pi1'].dataFromData( pi1 )
            
        for ic, (irc, icc) in enumerate( iterSym( dim ) ):
            omega2 = corrsRS.statesRS[irc,icc][indx]
            pi2 = pis[irc,icc] + omega2
            # Non-state variables must be assigned manually.
            problem.variables['Pi2'].dataFromData( pi2 )

            # Variables have their data, so evaluate the term.
            val = evalTermOP( None, coefTerm, problem )

            coef[ir,ic] = val
    return coef


##
# c: 05.05.2008, r: 05.05.2008
def main():
    from sfe.base.base import spause
    from sfe.base.conf import ProblemConf, getStandardKeywords
    from sfe.fem.problemDef import ProblemDefinition
    from sfe.base.ioutils import getTrunk

    nm.set_printoptions( precision = 3 )

    spause( r""">>>
First, this file will be read in place of an input
(problem description) file.
Press 'q' to quit the example, press any other key to continue...""" )
    required, other = getStandardKeywords()
    # Use this file as the input file.
    conf = ProblemConf.fromFile( __file__, required, other )
    print conf
    spause( r""">>>
...the read input.
['q'/other key to quit/continue...]""" )

    spause( r""">>>
Now the input will be used to create a ProblemDefinition instance.
['q'/other key to quit/continue...]""" )
    problem = ProblemDefinition.fromConf( conf,
                                          initVariables = False,
                                          initEquations = False )
    print problem
    spause( r""">>>
...the ProblemDefinition instance.
['q'/other key to quit/continue...]""" )


    spause( r""">>>
The homogenized elastic coefficient $E_{ijkl}$ is expressed
using $\Pi$ operators, computed now. In fact, those operators are permuted
coordinates of the mesh nodes.
['q'/other key to quit/continue...]""" )
    pis = createPis( problem, conf.variables, 'Pi' )
    print pis
    spause( r""">>>
...the $\Pi$ operators.
['q'/other key to quit/continue...]""" )

    ofnTrunk = getTrunk( conf.fileName_mesh ) + '_out'
    spause( r""">>>
Next, $E_{ijkl}$ needs so called steady state correctors $\bar{\omega}^{rs}$,
computed now. The results will be saved in: %s_*.vtk
['q'/other key to quit/continue...]""" % ofnTrunk )

    corrsRS = solveSteadyCorrectors_rs( problem, conf.equations,
                                        conf.variables, pis, ofnTrunk )
    print corrsRS
    spause( r""">>>
...the $\bar{\omega}^{rs}$ correctors.
['q'/other key to quit/continue...]""" )


    spause( r""">>>
Finally, $E_{ijkl}$ can be computed.
['q'/other key to quit/continue...]""" )
    cE = coefE( problem, corrsRS, pis )
    print r""">>>
The homogenized elastic coefficient $E_{ijkl}$, symmetric storage
with rows, columns in 11, 22, 12 ordering:"""
    print cE
    
if __name__ == '__main__':
    main()

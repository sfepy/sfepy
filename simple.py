#!/usr/bin/env python
# 12.01.2007, c 
import os.path as op
from optparse import OptionParser

import init_sfe
from sfe.base.base import *
from sfe.base.conf import ProblemConf
import sfe.base.ioutils as io
from sfe.fem.problemDef import ProblemDefinition
from sfe.solvers.generic import solveStationary, saveOnly

##
# c: 12.01.2007, r: 02.01.2008
def solveDirect( conf, options ):
    """Generic (simple) stationary problem solver."""
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
    
    dpb, vecDP, data = solveStationary( conf, saveNames = saveNames )
    out = dpb.stateToOutput( vecDP )

    fd = open( ofnTrunk + '.vtk', 'w' )
    io.writeVTK( fd, dpb.domain.mesh, out )
    fd.close()

    if options.dump:
        import tables as pt
        import numarray as nar

        fd = pt.openFile( ofnTrunk + '.h5', mode = 'w', title = "Dump file" )
        for key, val in out.iteritems():
            fd.createArray( fd.root, key, nar.asarray( val.data ), 
                            '%s data' % val.mode )
        fd.close()

    return dpb, vecDP, data

##
# 26.03.2007, c
def printTerms():
    import sfe.terms as t
    tt = t.termTable
    ct = t.cacheTable
    print 'Terms: %d available:' % len( tt )
    print sorted( tt.keys() )
    print 'Term caches: %d available:' % len( ct )
    print sorted( ct.keys() )

usage = """%prog [options] fileNameIn"""

help = {
    'fileName' :
    'basename of output file(s) [default: <basename of input file>]',
    'dump' :
    "dump problem state [default: %default]",
    'saveEBC' :
    "save problem state showing EBC (Dirichlet conditions) [default: %default]",
    'saveRegions' :
    "save problem regions as meshes [default: %default]",
    'saveFieldMeshes' :
    "save meshes of problem fields (with extra DOF nodes) [default: %default]",
    'saveRegionFieldMeshes' :
    "save meshes of regions of problem fields (with extra DOF nodes)"
    "[default: %default]",
    'solveNot' :
    "do not solve (use in connection with --save-*) [default: %default]",
    'list' :
    "list data according to what, what can be one of: terms",
}

##
# c: 12.01.2007, r: 23.01.2008
def main():
    version = open( op.join( init_sfe.install_dir,
                             'VERSION' ) ).readlines()[0][:-1]

    parser = OptionParser( usage = usage, version = "%prog " + version )
    parser.add_option( "-o", "", metavar = 'fileName',
                       action = "store", dest = "outputFileNameTrunk",
                       default = None, help = help['fileName'] )
    parser.add_option( "", "--dump",
                       action = "store_true", dest = "dump",
                       default = False, help = help['dump'] )
    parser.add_option( "", "--save-ebc",
                       action = "store_true", dest = "saveEBC",
                       default = False, help = help['saveEBC'] )
    parser.add_option( "", "--save-regions",
                       action = "store_true", dest = "saveRegions",
                       default = False, help = help['saveRegions'] )
    parser.add_option( "", "--save-field-meshes",
                       action = "store_true", dest = "saveFieldMeshes",
                       default = False, help = help['saveFieldMeshes'] )
    parser.add_option( "", "--save-region-field-meshes",
                       action = "store_true", dest = "saveRegionFieldMeshes",
                       default = False, help = help['saveRegionFieldMeshes'] )
    parser.add_option( "", "--solve-not",
                       action = "store_true", dest = "solveNot",
                       default = False, help = help['solveNot'] )
    parser.add_option( "", "--list", metavar = 'what',
                       action = "store", dest = "_list",
                       default = None, help = help['list'] )

    options, args = parser.parse_args()
#    print options; pause()

    if (len( args ) == 1):
        fileNameIn = args[0];
    else:
        if options._list == 'terms':
            printTerms()
        else:
            parser.print_help(),
        return
    
    required = ['fileName_mesh', 'field_[0-9]+', 'ebc|nbc', 'fe', 'equations',
                'region_[0-9]+|regions', 'variables', 'material_[0-9]+',
                'integral_[0-9]+', 'solver_[0-9]+']
    other = ['functions', 'modules', 'epbc', 'lcbc', 'options']
    if options.solveNot:
        required.remove( 'equations' )
        required.remove( 'solver_[0-9]+' )
        other += ['equations']

    conf = ProblemConf.fromFile( fileNameIn, required, other )
##     print conf
##     pause()

    dpb, vecDP, data = solveDirect( conf, options )

    opts = conf.options
    if hasattr( opts, 'postProcessHook' ): # User postprocessing.
        postProcessHook = getattr( conf, opts.postProcessHook )
        postProcessHook( dpb, vecDP, data )

if __name__ == '__main__':
    main()

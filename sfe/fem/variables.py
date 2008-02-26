from sfe.base.base import *
import sfe.base.la as la
from extmods.fem import rawGraph

isState = 0
isVirtual = 1
isParameter = 2
isOther = 3
isField = 10

##
# 11.07.2007, c
# 19.07.2007
def createDofConn( conn, dpn, imax ):
    if dpn == 1:
        dc = conn.copy()
    else:
        nEl, nEP = conn.shape
        nED = nEP * dpn
        dc = nm.empty( (nEl, nED), dtype = conn.dtype )
        for ic in range( nED ):
            inod = ic / dpn
            idof = ic % dpn
##                    iloc = ic
            iloc = nEP * idof + inod # Hack: For DBD order.
            dc[:,iloc] = dpn * conn[:,inod] + idof

    try:
        imax = max( imax, nm.amax( dc.ravel() ) )
    except: # empty dc (a non-first point dc - e.g. springs)
        pass

    return dc, imax

##
# c: 11.07.2007, r: 04.02.2008
def createADofConns( eq, iterator, indx ):
    adcs = {}
    for key, dc in iterator():
        if isinstance( dc, dict ):
            adcss = createADofConns( eq, dc.iteritems, indx )
            for subkey, subdc in adcss.iteritems():
                adcs[(key, subkey)] = subdc
        elif dc is not None:
            aux = eq[dc]
            adcs[key] = aux + nm.asarray( nm.where( aux >= 0, indx.start, 0 ),
					  dtype = nm.int32 )
        else:
            adcs[key] = None
    return adcs

##
# 26.07.2006, c
def zeroConfEBC( conf ):
    new = {}
    for rname, bcs in conf.iteritems():
        newbc = []
        for bc in bcs:
            newbc.append( bc[:-1] + (0.0,) )
        new[rname] = tuple( newbc )
    return new

##
# 27.11.2006, c
def invertDataShapes( varShapes ):
    shapes = {}
    for name, groups in varShapes.iteritems():
        for ig, shape in groups.iteritems():
            shapes.setdefault( ig, {} )[name] = shape
##                 if not shapes.has_key( ig ):
##                     shapes[ig] = {}
##                 shapes[ig][name] = shape
    return shapes

##
# 15.03.2007, c
# 04.06.2007
def resolveChains( masterSlave, chains ):
    """Resolve EPBC chains - e.g. in corner nodes."""

    for chain in chains:
        slave = chain[-1]
        masterSlave[chain[:-1]] = slave + 1
        masterSlave[slave] = - chain[0] - 1 # Any of masters...

##
# 04.06.2007, c
def groupChains( chainList ):

    chains = []
    while len( chainList ):
        chain = set( chainList.pop( 0 ) )
##         print ':', chain
        ii = 0
        while ii < len( chainList ):
            c1 = sorted( chainList[ii] )
#            print '--', ii, c1, chain
            is0 = c1[0] in chain
            is1 = c1[1] in chain

            if is0 and is1:
                chainList.pop( ii )
            elif is0 or is1:
                chain.update( c1 )
                chainList.pop( ii )
                ii = 0
            else:
                ii += 1
#            print ii, chain, chainList
##         print '->', chain
##         print chainList
##         pause()

        chains.append( list( chain ) )

#    print 'EPBC chain groups:', chains
    aux = {}
    for chain in chains:
        aux.setdefault( len( chain ), [0] )[0] += 1
#    print 'EPBC chain counts:', aux

    return chains

##
# 19.07.2006
class DofInfo( Struct ):
    pass

##
# 14.07.2006, c
class Variables( Container ):

    ##
    # c: 14.07.2006, r: 18.02.2008
    def fromConf( conf, fields ):
        objs = OneTypeList( Variable )
        state, virtual, parameter, other = [], [], [], []
        for ii, (key, val) in enumerate( conf.iteritems() ):
            var = Variable.fromConf( key, val, fields )
            if var.isState():
                state.append( ii )
            elif var.isVirtual():
                virtual.append( ii )
            elif var.isParameter():
                parameter.append( ii )
            elif var.isOther():
                other.append( ii )
            objs.append( var )

        obj = Variables( objs,
                         state = state,
                         virtual = virtual,
                         parameter = parameter,
                         other = other,
                         domain = fields[0].domain,
                         fields = fields,
                         hasVirtualDCs = False,
                         hasLCBC = False )

        indx = nm.array( [var._order for var in obj.iterState()] )
        obj.orderedState = indx.argsort()

        obj.linkDuals()
        return obj
    fromConf = staticmethod( fromConf )

    ##
    # 05.09.2007, c
    def linkDuals( self ):
        for ii in self.virtual:
            vvar = self[ii]
            self[vvar.primaryVarName].dualVarName = vvar.name

    ##
    # 26.07.2007, c
    def getNames( self, kind = None ):
        if kind is None:
            names = [var.name for var in self]
        else:
            names = [var.name for var in self if var.isKind( kind )]
        return names

    ##
    # 07.10.2005, c
    # 25.10.2005
    # 26.10.2005
    # 01.11.2005
    # 04.11.2005
    # 07.03.2006
    # 19.07.2006
    # 12.07.2007
    # 05.09.2007
    def setupDofInfo( self, makeVirtual = False ):
        """Sets also iDofMax."""
        def _setupDofInfo( iterable ):
            ptr = [0]
            nNod = []
            dpn = []
            vnames = []
            for ii in iterable:
                var = self[ii]
                var.iDofMax = var.field.nNod * var.dpn

                nNod.append( var.field.nNod )
                dpn.append( var.dpn )
                ptr.append( ptr[-1] + dpn[-1] * nNod[-1] )
                vnames.append( var.name )

            di = DofInfo(
                name = 'dofInfo',
                ptr = nm.array( ptr, dtype = nm.int32 ),
                nNod = nm.array( nNod, dtype = nm.int32 ),
                dpn = nm.array( dpn, dtype = nm.int32 ),
                vnames = vnames, indx = {}, nDofs = {}
            )

            for ii, name in enumerate( di.vnames ):
                di.indx[name] = slice( di.ptr[ii], di.ptr[ii+1] )
                di.nDofs[name] = di.ptr[ii+1] - di.ptr[ii]
            return di

        self.di = _setupDofInfo( self.state )
        if makeVirtual:
            self.vdi = _setupDofInfo( self.virtual )
        else:
            self.vdi = self.di

    ##
    # c: 16.10.2006, r: 25.02.2008
    def _listBCOfVars( self, bcDefs, isEBC = True ):

        bcOfVars = dictFromKeysInit( (key for key in self.di.vnames), list )
        if bcDefs is None: return bcOfVars

        for key, bc in bcDefs.iteritems():
##             print key
##             print bc
            for dofs, val in bc.dofs.iteritems():
                vname = dofs.split( '.' )[0]
                vbc = copy( bc )
                vbc.dofs = (dofs, val)
##                 print '///', vbc
                bcOfVars[vname].append( (key, vbc) )

        return bcOfVars

    ##
    # c: 03.10.2007, r: 18.02.2008
    def setupLCBCOperators( self, lcbc, regions ):
        if lcbc is None: return

        self.hasLCBC =True

        lcbcOfVars = self._listBCOfVars( lcbc )

        # Assume disjoint regions.
        lcbcOps = {}
        for varName, bcs in lcbcOfVars.iteritems():
            var = self[varName]
            lcbcOps[varName] = var.createLCBCOperators( bcs, regions )

        opsLC = []
        eqLCBC = nm.empty( (0,), dtype = nm.float64 )
        nGroups = 0
        for varName, lcbcOp in lcbcOps.iteritems():
            if lcbcOp is None: continue
#            print varName, lcbcOp

            indx = self.adi.indx[varName]
            aux = nm.where( lcbcOp.eqLCBC >= 0, indx.start, 0 )
            eqLCBC = nm.hstack( (eqLCBC, lcbcOp.eqLCBC + aux) )
            opsLC.extend( lcbcOp.opsLC )

            nRigidDof = lcbcOp.nRigidDof
            dim = lcbcOp.dim
            nGroups += lcbcOp.nGroups

        if nGroups == 0:
            self.hasLCBC = False
            return
            
        nDof = self.adi.ptr[-1]

        ii = nm.nonzero( eqLCBC )[0]
        nConstrained = ii.shape[0]
        nDofNotRigid = nDof - nConstrained
        nDofReduced = nDofNotRigid + nGroups * nRigidDof
        print nDof, nDofReduced, nConstrained, nDofNotRigid

        mtxLC = sp.lil_matrix( (nDof, nDofReduced), dtype = nm.float64 )
        ir = nm.where( eqLCBC == 0 )[0]
        ic = nm.arange( nDofReduced, dtype = nm.int32 )
        mtxLC[ir,ic] = 1.0
        for ii, opLC in enumerate( opsLC ):
            indx = nm.where( eqLCBC == (ii + 1) )[0]
            icols = slice( nDofNotRigid + nRigidDof * ii,
                           nDofNotRigid + nRigidDof * (ii + 1) )
            mtxLC[indx,icols] = opLC

        mtxLC = mtxLC.tocsr()
##         import pylab
##         from sfe.base.plotutils import spy
##         spy( mtxLC )
##         pylab.show()
##         print mtxLC
        nnz = nDof - nConstrained + nConstrained * dim
        print nnz, mtxLC.getnnz()
        assert nnz >= mtxLC.getnnz()

        self.opLCBC = mtxLC

    ##
    # 04.10.2007, c
    def getLCBCOperator( self ):
        if self.hasLCBC:
            return self.opLCBC
        else:
            print 'no LCBC defined!'
            raise ValueError

    ##
    # c: 01.11.2005, r: 18.02.2008
    def equationMapping( self, ebc, epbc, regions, ts, funmod,
                         vregions = None ):

        if vregions is None:
            vregions = regions

        ##
        # Assing EBC, PBC to variables and regions.
        self.bcOfVars = self._listBCOfVars( ebc )
        dictExtend( self.bcOfVars, self._listBCOfVars( epbc, isEBC = False ) )

        ##
        # List EBC nodes/dofs for each variable.
        for varName, bcs in self.bcOfVars.iteritems():
            var = self[varName]
            var.equationMapping( bcs, regions, self.di, ts, funmod )
            if self.hasVirtualDCs:
                vvar = self[var.dualVarName]
                vvar.equationMapping( bcs, vregions, self.vdi, ts, funmod )

##             print var.eqMap
##             pause()

        ##
        # Adjust by offsets - create active dof info.
        def _createADofInfo( di ):
            adi = DofInfo(
                name = 'active_dofInfo',
                ptr = nm.array( di.ptr, dtype = nm.int32 ),
                nNod = nm.array( di.nNod, dtype = nm.int32 ),
                dpn = nm.array( di.dpn, dtype = nm.int32 ),
                vnames = di.vnames, indx = {}, nDofs = {}
            )
            for ii, key in enumerate( adi.vnames ):
                adi.nDofs[key] = self[key].eqMap.nEq
                adi.ptr[ii+1] = adi.ptr[ii] + adi.nDofs[key]
                adi.indx[key] = slice( adi.ptr[ii], adi.ptr[ii+1] )
            return adi

        self.adi = _createADofInfo( self.di )
        if self.hasVirtualDCs:
            self.avdi = _createADofInfo( self.vdi )
        else:
            self.avdi = self.adi

    ##
    # c: 09.01.2008, r: 09.01.2008
    def getNodesOfGlobalDofs( self, igdofs ):
        """not stripped..."""
        di = self.di
        
        nods = nm.empty( (0,), dtype = nm.int32 )
        for ii in self.state:
            var = self[ii]
            indx = di.indx[var.name]
            igdof = igdofs[(igdofs >= indx.start) & (igdofs < indx.stop)]
            ivdof = igdof - indx.start
            inod = ivdof / var.dpn
            nods = nm.concatenate( (nods, inod) )
##             print var.name, indx
##             print igdof
##             print ivdof
##             print inod
##             pause()
        return nods

    ##
    # c: 23.11.2005, r: 26.02.2008
    def setupDofConns( self, makeVirtual = False ):
        output( 'setting up dof connectivities...' )
        tt = time.clock()

        for ii in self.state:
            var = self[ii]
            var.setupDofConns()

        if makeVirtual:
            for ii in self.virtual:
                var = self[ii]
                var.setupDofConns()
            self.hasVirtualDCs = True
        else:
            self.hasVirtualDCs = False

        output( '...done in %.2f s' % (time.clock() - tt) )

    ##
    # 08.08.2006, c
    # 11.10.2006
    # 20.02.2007
    # 22.02.2007
    # 11.07.2007
    # 05.09.2007
    def setupADofConns( self ):
        """Translate dofs to active dofs."""
        def _setupADofConns( iterable, adi ):
            adofConns = dictFromKeysInit( [self[ii].name for ii in iterable],
                                          Struct )
            for ii in iterable:
                var = self[ii]
                indx = adi.indx[var.name]
                eq = var.eqMap.eq
                adofConns[var.name].name = 'adofConns'
                it =  var.iterDofConns( 'volume' )
                adofConns[var.name].volumeDCs = createADofConns( eq, it, indx )
                it =  var.iterDofConns( 'surface' )
                adofConns[var.name].surfaceDCs = createADofConns( eq, it, indx )
                it =  var.iterDofConns( 'edge' )
                adofConns[var.name].edgeDCs = createADofConns( eq, it, indx )
                it =  var.iterDofConns( 'point' )
                adofConns[var.name].pointDCs = createADofConns( eq, it, indx )
            return adofConns

        self.adofConns = _setupADofConns( self.state, self.adi )
        if self.hasVirtualDCs:
            self.avdofConns = _setupADofConns( self.virtual, self.avdi )
        else:
            self.avdofConns = self.adofConns

##         print self.adofConns.values()[0]
##         pause()

    ##
    # c: 10.12.2007, r: 15.01.2008
    def getADofConn( self, varName, isDual, dcType, ig ):
        """Note that primary and dual variables must have same Region!"""
        kind, regionName = dcType

        var = self[varName]
        if isDual:
            if not self.hasVirtualDCs:
                varName = var.primaryVarName
            adcs = self.avdofConns[varName]
        else:
            adcs = self.adofConns[varName]

        if kind == 'volume':
            try:
                dc = adcs.volumeDCs[ig]
            except:
                debug()
        else:
            if kind == 'surface':
                dcs = adcs.surfaceDCs
            elif kind == 'edge':
                dcs = adcs.edgeDCs
            elif kind == 'point':
                dcs = adcs.pointDCs
            else:
                print 'uknown dof connectivity kind:', kind
                raise ValueError
            dc = dcs[(ig, regionName)]
        return dc
        
    ##
    # 05.09.2007, c
    def fixCoarseGridADofConns( self, iemaps, varName ):
        """Volume only!"""
        dcs = self.adofConns[varName].volumeDCs
        for ig, dc in dcs.iteritems():
            iemap = iemaps[ig]
            dcs[ig] = dc[iemap]

    ##
    # c: 23.11.2005, r: 04.02.2008
    def createMatrixGraph( self, varNames = None, vvarNames = None ):
        """
        Create tangent matrix graph. Order of dof connectivities is not
        important here...
        """
        def _prepareDCLists( adofConns, varNames = None ):
            if varNames is None:
                varNames = adofConns.iterkeys()

            gdcs = {}
            for varName in varNames:
                adcs = adofConns[varName]
                for ig, dc in adcs.volumeDCs.iteritems():
##                     print dc
                    gdcs.setdefault( ig, [] ).append( dc )
            return gdcs

        shape = (self.avdi.ptr[-1], self.adi.ptr[-1])
        output( 'matrix shape:', shape )

        cgdcs = _prepareDCLists( self.adofConns, varNames )
##         print cgdcs
##         pause()
        if self.hasVirtualDCs:
            rgdcs = _prepareDCLists( self.avdofConns, vvarNames )
        else:
            rgdcs = cgdcs

        ##
        # Make all permutations per element group.
        rdcs = []
        cdcs = []
        for ig in rgdcs.iterkeys():
            rgdc, cgdc = rgdcs[ig], cgdcs[ig]
            for perm in la.cycle( [len( rgdc ), len( cgdc )] ):
                rdcs.append( rgdc[perm[0]] )
                cdcs.append( cgdc[perm[1]] )
#                print ' ', perm, '->', rdcs[-1].shape, cdcs[-1].shape

        output( 'assembling matrix graph...' )
        tt = time.clock()

#	shape = nm.array( shape, dtype = nm.long )
        ret, prow, icol = rawGraph( int( shape[0] ), int( shape[1] ),
                                    len( rdcs ), rdcs, cdcs )
        output( '...done in %.2f s' % (time.clock() - tt) )
        nnz = prow[-1]
##         print ret, prow, icol, nnz
	
        data = nm.zeros( (nnz,), dtype = nm.float64 )
        matrix = sp.csr_matrix( (data, icol, prow), shape )
##         matrix.save( 'matrix', format = '%d %d %e\n' )
##         pause()

        return matrix

    ##
    # 24.07.2006, c
    # 25.07.2006
    # 04.08.2006
    def dataFromState( self, state = None ):
        for ii in self.state:
            var = self[ii]
            var.dataFromState( state, self.di.indx[var.name] )

    ##
    # 25.07.2006, c
    def createStateVector( self ):
        vec = nm.zeros( (self.di.ptr[-1],), nm.float64 )
        return vec

    ##
    # 25.07.2006, c
    def createStrippedStateVector( self ):
        vec = nm.zeros( (self.adi.ptr[-1],), nm.float64 )
        return vec

    ##
    # 22.11.2005, c
    # 25.07.2006
    # 19.09.2006
    # 18.10.2006
    def applyEBC( self, vec, forceValues = None ):
        for varName in self.bcOfVars.iterkeys():
            eqMap = self[varName].eqMap
            i0 = self.di.indx[varName].start
            ii = i0 + eqMap.eqEBC
##             print ii, eqMap.valEBC
##             pause()
            if forceValues is None:
                vec[ii] = eqMap.valEBC
            else:
                if isinstance( forceValues, dict ):
                    vec[ii] = forceValues[varName]
                else:
                    vec[ii] = forceValues
            # EPBC.
            vec[i0+eqMap.master] = vec[i0+eqMap.slave]

    ##
    # 27.11.2005, c
    # 09.12.2005
    # 25.07.2006
    # 18.10.2006
    def updateVec( self, vec, delta ):
        for varName in self.bcOfVars.iterkeys():
            eqMap = self[varName].eqMap
            i0 = self.di.indx[varName].start
            ii = i0 + eqMap.eqi
##            print ii.shape, delta[adi.indx[varName]].shape
            vec[ii] -= delta[self.adi.indx[varName]]
            # EPBC.
            vec[i0+eqMap.master] = vec[i0+eqMap.slave]

    ##
    # 12.04.2007, c
    def stripStateVector( self, vec ):
        svec = nm.empty( (self.adi.ptr[-1],), nm.float64 )
        for varName in self.bcOfVars.iterkeys():
            eqMap = self[varName].eqMap
            i0 = self.di.indx[varName].start
            ii = i0 + eqMap.eqi
##            print ii.shape, delta[adi.indx[varName]].shape
            svec[self.adi.indx[varName]] = vec[ii]
        return svec

    ##
    # created:       12.04.2007
    # last revision: 14.12.2007
    def makeFullVec( self, svec, varName = None, forceValue = None ):
        """
        Make a full vector satisfying E(P)BC
        from a stripped vector. For a selected variable if varName is set.
        """
        def _makeFullVec( vec, eqMap ):
            # EBC.
            ii = eqMap.eqEBC
            if forceValue is None:
                vec[ii] = eqMap.valEBC
            else:
                vec[ii] = forceValue

            # Stripped vector values.
            ii = eqMap.eqi
            vec[ii] = svec

            # EPBC.
            vec[eqMap.master] = vec[eqMap.slave]

        if self.hasLCBC:
            svec = self.opLCBC * svec

        if varName is None:
            vec = self.createStateVector()

            for varName in self.bcOfVars.iterkeys():
                eqMap = self[varName].eqMap
                _makeFullVec( vec[self.di.indx[varName]], eqMap )
        else:
            vec = nm.empty( (self.di.nDofs[varName],), nm.float64 )
            eqMap = self[varName].eqMap
            _makeFullVec( vec, eqMap )

        return vec

    ##
    # 14.03.2007, c
    def hasEBC( self, vec, forceValues = None ):
        for varName in self.bcOfVars.iterkeys():
            eqMap = self[varName].eqMap
            i0 = self.di.indx[varName].start
            ii = i0 + eqMap.eqEBC
            if forceValues is None:
                if not nm.allclose( vec[ii], eqMap.valEBC ):
                    return False
            else:
                if isinstance( forceValues, dict ):
                    if not nm.allclose( vec[ii], forceValues[varName] ):
                        return False
                else:
                    if not nm.allclose( vec[ii], forceValues ):
                        return False
            # EPBC.
            if not nm.allclose( vec[i0+eqMap.master], vec[i0+eqMap.slave] ):
                return False
        return True

    ##
    # 26.07.2007, c
    def getIndx( self, varName, stripped = False ):
        if self[varName].isState():
            if stripped:
                return self.adi.indx[varName]
            else:
                return self.di.indx[varName]
        else:
            output( '%s is not a state part' % varName )
            raise IndexError

    ##
    # 26.07.2006, c
    # 12.04.2007
    # 26.07.2007
    def getStatePartView( self, state, varName, stripped = False ):
        return state[self.getIndx( varName, stripped )]

    ##
    # 26.07.2006, c
    # 12.04.2007
    # 26.07.2007
    def setStatePart( self, state, part, varName, stripped = False ):
        state[self.getIndx( varName, stripped )] = part

    ##
    # 26.07.2006, c
    # 02.08.2006
    # 04.08.2006
    def nonStateDataFromState( self, varNames, state, varNamesState ):
        if isinstance( varNames, str ):
            varNames = [varNames]
            varNamesState = [varNamesState]

        for ii, varName in enumerate( varNames ):
            varNameState = varNamesState[ii]
            if self[varNameState].isState():
                self[varName].dataFromData( state,
                                            self.di.indx[varNameState] )
            else:
                output( '%s is not a state part' % varNameState )
                raise IndexError


    ##
    # Works for vertex data only.
    # c: 15.12.2005, r: 25.02.2008
    def stateToOutput( self, vec, fillValue = None ):

        nNod, di = self.domain.shape.nNod, self.di

        out = {}
        for key, indx in di.indx.iteritems():
            dpn = di.dpn[di.vnames.index( key )]
            aux = nm.reshape( vec[indx], (di.nDofs[key] / dpn, dpn) )
#            print key, aux.shape
            ext = self[key].extendData( aux, nNod, fillValue )
#            print ext.shape
#            pause()
            out[key] = Struct( name = 'output_data',
                               mode = 'vertex', data = ext,
                               dofs = self[key].dofs )
        return out

    ##
    # 24.08.2006, c
    # 20.09.2006
    def getFieldsOfVars( self, varNames ):
        fieldNames = {}
        for varName in varNames:
            if not self.has_key( varName ):
                raise RuntimeError, 'undefined variable %s' % varName
            fieldNames[varName] = self[varName].field.name
        return fieldNames

    ##
    # 27.11.2006, c
    def iterState( self, ordered = False ):

        if ordered:
            for ii in self.orderedState:
                yield self[self.state[ii]]

        else:
            for ii in self.state:
                yield self[ii]

##
# 11.07.2006, c
class Variable( Struct ):

    ##
    # c: 14.07.2006?, r: 25.02.2008
    def fromConf( key, conf, fields ):
        flags = set()
        kind, family = conf.kind.split()

        obj = Variable( flags, name = conf.name, key = key,
                        kind = kind, family = family )

        if kind == 'unknown':
            obj.flags.add( isState )
            if hasattr( conf, 'order' ):
                obj._order = int( conf.order )
            else:
                output( 'unnown variable %s: order missing' % conf.name )
                raise ValueError
            obj.dofName = obj.name
        elif kind == 'test':
            obj.flags.add( isVirtual )
            if hasattr( conf, 'dual' ):
                obj.primaryVarName = conf.dual
            else:
                output( 'test variable %s: related unknown missing' % conf.name )
                raise ValueError
            obj.dofName = obj.primaryVarName
        elif kind == 'parameter':
            obj.flags.add( isParameter )
            if hasattr( conf, 'like' ):
                obj.primaryVarName = conf.like
            else:
                output( 'parameter variable %s: related unknown missing'\
                        % conf.name )
                raise ValueError
            obj.dofName = obj.primaryVarName
        else:
            obj.flags.add( isOther )
            print 'unknown variable family: %s' % family
            raise NotImplementedError

        if family == 'field':
            try:
                fld = fields[conf.field]
            except:
                output( 'field "%s" does not exist!' % conf.fields )
                raise

            obj.setField( fld )

        return obj
    fromConf = staticmethod( fromConf )

    ##
    # 11.07.2006, c
    # 04.08.2006
    def __init__( self, flags, data = None, indx = 0, **kwargs ):
        Struct.__init__( self, **kwargs )

        self.flags = set()
        for flag in flags:
            self.flags.add( flag )

        self.data = data
        self.indx = None
        self.nDof = None
        self.currentAp = None

        if self.isVirtual():
            self.data = None

    ##
    # 11.07.2006, c
    def __call__( self ):
        return (self.data, self.indx)

    ##
    # 11.07.2006, c
    def isState( self ):
        return isState in self.flags

    ##
    # 11.07.2006, c
    def isVirtual( self ):
        return isVirtual in self.flags

    ##
    # 26.07.2007, c
    def isParameter( self ):
        return isParameter in self.flags
    ##
    # 26.07.2007, c
    def isOther( self ):
        return isOther in self.flags

    ##
    # 26.07.2007, c
    def isKind( self, kind ):
        return eval( 'self.is%s()' % kind.capitalize() )

    ##
    # 26.07.2006, c
    def isNonStateField( self ):
        return (isField in self.flags)\
               and not (self.isState() or self.isVirtual())

    ##
    # c: 11.07.2006, r: 04.02.2008
    def dataFromState( self, state = None, indx = None ):
        if (not self.isState()) or (state is None): return

        self.data = state
        self.indx = slice( int( indx.start ), int( indx.stop ) )
        self.nDof = indx.stop - indx.start

    ##
    # c: 26.07.2006, r: 13.02.2008
    def dataFromData( self, data = None, indx = None ):
        if (not self.isNonStateField()) or (data is None): return

        self.data = data
        if indx is None:
            self.indx = slice( 0, len( data ) )
        else:
            self.indx = slice( int( indx.start ), int( indx.stop ) )
        self.nDof = self.indx.stop - self.indx.start

    ##
    # c: 11.07.2006, r: 25.02.2008
    def setField( self, field ):
        """Takes reference to a Field instance."""
        self.field = field
        self.dpn = nm.product( field.dim )

        if self.dofName is None:
            dofName = 'aux'
        else:
            dofName = self.dofName
        self.dofs = [dofName + ('.%d' % ii) for ii in range( self.dpn )]

        self.flags.add( isField )

    ##
    # c: 08.08.2006, r: 15.01.2008
    def setupDofConns( self ):
        dpn = self.dpn
        field = self.field

        dofConns = Struct( name = 'dofConns', volumeDCs = {},
                           surfaceDCs = {}, edgeDCs = {}, pointDCs = {} )
        imax = -1
        ##
        # Expand nodes into dofs.
        for regionName, ig, ap in field.aps.iterAps():

            dc, imax = createDofConn( ap.econn, dpn, imax )
            dofConns.volumeDCs[ig] = dc

            if ap.surfaceData:
                dcs = {}
                for key, sd in ap.surfaceData.iteritems():
                    dc, imax2 = createDofConn( sd.econn, dpn, 0 )
                    assert imax2 <= imax
                    dcs[key] = dc
                dofConns.surfaceDCs[ig] = dcs

            else:
                dofConns.surfaceDCs[ig] = None

            if ap.edgeData:
                raise NotImplementedError
            else:
                dofConns.edgeDCs[ig] = None

            if ap.pointData:
                dcs = {}
                for key, conn in ap.pointData.iteritems():
                    dc, imax2 = createDofConn( conn, dpn, 0 )
                    # imax2 can be greater than imax, as all spring nodes
                    # are assigned to the first group!!!
##                     print conn
##                     print dc
##                     print key, dc.shape
##                     pause()
                    dcs[key] = dc
                dofConns.pointDCs[ig] = dcs

            else:
                dofConns.pointDCs[ig] = None

##         print dofConns
##         pause()
        self.dofConns = dofConns
        assert self.iDofMax >= imax

    ##
    # 20.02.2007, c
    # 11.07.2007
    def iterDofConns( self, kind ):
        if kind == 'volume':
            return self.dofConns.volumeDCs.iteritems
        elif kind == 'surface':
            return self.dofConns.surfaceDCs.iteritems
        elif kind == 'edge':
            return self.dofConns.edgeDCs.iteritems
        elif kind == 'point':
            return self.dofConns.pointDCs.iteritems
        else:
            print 'uknown dof connectivity kind:', kind
            raise ValueError

    ##
    # c: 25.02.2008, r: 25.02.2008
    def _canonize( self, dofs ):
        vname, dd = dofs.split( '.' )
        if dd == 'all':
            cdofs = self.dofs
        elif dd[0] == '[':
            cdofs = [vname + '.' + ii.strip()
                     for ii in dd[1:-1].split( ',' )]
        else:
            cdofs = [dofs]
        return cdofs

    ##
    # c: 18.10.2006, r: 25.02.2008
    def expandNodesToEquations( self, nods, dofs ):

        eq = nm.array( [], dtype = nm.int32 )
        for dof in self._canonize( dofs ):
            idof = self.dofs.index( dof )
            eq = nm.concatenate( (eq, self.dpn * nods + idof) )
        return eq

    ##
    # c: 03.10.2007, r: 25.02.2008
    def createLCBCOperators( self, bcs, regions ):
        if len( bcs ) == 0: return None

        eq = self.eqMap.eq
        nDof = self.eqMap.nEq
        
        nGroups = len( bcs )
        eqLCBC = nm.zeros( (nDof,), dtype = nm.int32 )
        
        opsLC = []
        for ii, (key, bc) in enumerate( bcs ):
            print self.name, bc.name

            region = regions[bc.region]
            print region.name

            nmaster = region.getFieldNodes( self.field, merge = True )
            print nmaster.shape

            dofs, kind = bc.dofs
            meq = eq[self.expandNodesToEquations( nmaster, dofs )]
            assert nm.all( meq >= 0 )
            
            eqLCBC[meq] = ii + 1
##             print meq, meq.shape
##             print nm.where( eqLCBC )[0]
            
            mcoor = self.field.getCoor( nmaster )[:,:-1]
            nNod, dim = mcoor.shape

#            print mcoor, mcoor.shape

            mtxE = nm.tile( nm.eye( dim, dtype = nm.float64 ), (nNod, 1) )
            if dim == 2:
                mtxR = nm.empty( (dim * nNod, 1), dtype = nm.float64 )
                mtxR[0::dim,0] = -mcoor[:,1]
                mtxR[1::dim,0] = mcoor[:,0]
                nRigidDof = 3
            elif dim == 3:
                mtxR = nm.zeros( (dim * nNod, dim), dtype = nm.float64 )
                mtxR[0::dim,1] = mcoor[:,2]
                mtxR[0::dim,2] = -mcoor[:,1]
                mtxR[1::dim,0] = -mcoor[:,2]
                mtxR[1::dim,2] = mcoor[:,0]
                mtxR[2::dim,0] = mcoor[:,1]
                mtxR[2::dim,1] = -mcoor[:,0]
                nRigidDof = 6
            else:
                print 'dimension in [2,3]: %d' % dim
                raise ValueError

            opLC = nm.hstack( (mtxR, mtxE) )
##             print opLC, opLC.shape

            opsLC.append( opLC )

        return Struct( eqLCBC = eqLCBC,
                       opsLC = opsLC,
                       nGroups = nGroups,
                       nRigidDof = nRigidDof,
                       dim = dim )

    ##
    # c: 20.07.2006, split, r: 25.02.2008
    def equationMapping( self, bcs, regions, di, ts, funmod, warn = False ):
        """EPBC: master and slave dofs must belong to the same field (variables
        can differ, though)."""
        # Sort by ebc definition name.
        bcs.sort( cmp = lambda i1, i2: cmp( i1[0], i2[0] ) )
##         print bcs
##         pause()
        
        self.eqMap = eqMap = Struct()

        eqMap.eq = nm.arange( di.nDofs[self.name], dtype = nm.int32 )
        eqMap.valEBC = nm.empty( (0,), dtype = nm.float64 )
        if len( bcs ) == 0:
            ##
            # No ebc for this field.
            eqMap.eqi = nm.arange( di.nDofs[self.name], dtype = nm.int32 )
            eqMap.eqEBC = nm.empty( (0,), dtype = nm.int32 )
            eqMap.nEq = eqMap.eqi.shape[0]
            eqMap.nEBC = eqMap.eqEBC.shape[0]
            eqMap.master = nm.empty( (0,), dtype = nm.int32 )
            eqMap.slave = nm.empty( (0,), dtype = nm.int32 )
            return

        field = self.field

        eqEBC = nm.zeros( (di.nDofs[self.name],), dtype = nm.int32 )
        valEBC = nm.zeros( (di.nDofs[self.name],), dtype = nm.float64 )
        masterSlave = nm.zeros( (di.nDofs[self.name],), dtype = nm.int32 )
        chains = []
        for key, bc in bcs:
            if key[:3] == 'ebc':
                ntype = 'EBC'
                rname = bc.region
            else:
                ntype = 'EPBC'
                rname = bc.region[0]

            try:
                region = regions[rname]
            except KeyError:
                print "no region '%s' used in BC %s!" % (rname, bc)
                raise

##             print ir, key, bc
##             debug()
            # Get master region nodes.
            masterNodList = region.getFieldNodes( field )
            for master in masterNodList[:]:
                if master is None:
                    masterNodList.remove( master )
                    if warn:
                        output( 'warning: ignoring nonexistent %s'\
                                + ' node (%s) in %s'\
                                % (ntype, self.name, region.name) )

            if len( masterNodList ) == 0:
                continue

            if ntype == 'EBC': # EBC.
                dofs, val = bc.dofs
                ##
                # Evaluate EBC values.
                vv = nm.empty( (0,), dtype = nm.float64 )
                nods = nm.unique1d( nm.hstack( masterNodList ) )
                coor = field.getCoor( nods )
                if type( val ) == str:
                    fun = getattr( funmod, val )
                    vv = fun( bc, ts, coor )
                else:
                    vv = nm.repeat( [val], nods.shape[0] * len( dofs ) )
##                 print nods
##                 print vv
                eq = self.expandNodesToEquations( nods, dofs )

                # Duplicates removed here...
                eqEBC[eq] = 1
                if vv is not None: valEBC[eq] = vv

            else: # EPBC.
                region = regions[bc.region[1]]
                slaveNodList = region.getFieldNodes( field )
                for slave in slaveNodList[:]:
                    if slave is None:
                        slaveNodList.remove( slave )
                        if warn:
                            output( 'warning: ignoring nonexistent EPBC'\
                                    + ' slave node (%s) in %s'\
                                    % (self.name, region.name) )
                if len( slaveNodList ) == 0:
                    continue

##                 print masterNodList
##                 print slaveNodList

                nmaster = nm.unique1d( nm.hstack( masterNodList ) )
                nslave = nm.unique1d( nm.hstack( slaveNodList ) )
##                 print nmaster + 1
##                 print nslave + 1
                if nmaster.shape != nslave.shape:
                    raise 'EPBC list lengths do not match!\n(%s,\n %s)' %\
                          (nmaster, nslave)

                mcoor = field.getCoor( nmaster )
                scoor = field.getCoor( nslave )
                fun = getattr( funmod, bc.match )
                i1, i2 = fun( mcoor, scoor )
##                 print nm.c_[mcoor[i1], scoor[i2]]
##                 print nm.c_[nmaster[i1], nslave[i2]] + 1

                meq = self.expandNodesToEquations( nmaster[i1], bc.dofs[0] )
                seq = self.expandNodesToEquations( nslave[i2], bc.dofs[1] )

                m_assigned = nm.where( masterSlave[meq] != 0 )[0]
                s_assigned = nm.where( masterSlave[seq] != 0 )[0]
                if m_assigned.size or s_assigned.size: # Chain EPBC.
##                     print m_assigned, meq[m_assigned]
##                     print s_assigned, seq[s_assigned]

                    aux = masterSlave[meq[m_assigned]]
                    sgn = nm.sign( aux )
                    om_chain = zip( meq[m_assigned], (aux - sgn) * sgn )
##                     print om_chain
                    chains.extend( om_chain )

                    aux = masterSlave[seq[s_assigned]]
                    sgn = nm.sign( aux )
                    os_chain = zip( seq[s_assigned], (aux - sgn) * sgn )
##                     print os_chain
                    chains.extend( os_chain )

                    m_chain = zip( meq[m_assigned], seq[m_assigned] )
##                     print m_chain
                    chains.extend( m_chain )

                    msd = nm.setdiff1d( s_assigned, m_assigned )
                    s_chain = zip( meq[msd], seq[msd] )
##                     print s_chain
                    chains.extend( s_chain )

                    msa = nm.union1d( m_assigned, s_assigned )
                    ii = nm.setdiff1d( nm.arange( meq.size ), msa )
                    masterSlave[meq[ii]] = seq[ii] + 1
                    masterSlave[seq[ii]] = - meq[ii] - 1
                else:
                    masterSlave[meq] = seq + 1
                    masterSlave[seq] = - meq - 1
##                 print 'ms', masterSlave
##                 print chains

##         print masterSlave
        chains = groupChains( chains )
        resolveChains( masterSlave, chains )

        ii = nm.argwhere( eqEBC == 1 )
        eqMap.eqEBC = nm.atleast_1d( ii.squeeze() )
        eqMap.valEBC = nm.atleast_1d( valEBC[ii].squeeze() )
        eqMap.master = nm.argwhere( masterSlave > 0 ).squeeze()
        eqMap.slave = masterSlave[eqMap.master] - 1
##         print eqMap.master
##         print eqMap.slave
##         pause()

        assert (eqMap.eqEBC.shape == eqMap.valEBC.shape)
##         print eqMap.eqEBC.shape
##         pause()
        eqMap.eq[eqMap.eqEBC] = -2
        eqMap.eq[eqMap.master] = -1
        eqMap.eqi = nm.compress( eqMap.eq >= 0, eqMap.eq )
        eqMap.eq[eqMap.eqi] = nm.arange( eqMap.eqi.shape[0], dtype = nm.int32 )
        eqMap.eq[eqMap.master] = eqMap.eq[eqMap.slave]
        eqMap.nEq = eqMap.eqi.shape[0]
        eqMap.nEBC = eqMap.eqEBC.shape[0]
        eqMap.nEPBC = eqMap.master.shape[0]
##         print eqMap
##         pause()

    ##
    # c: 24.07.2006, r: 15.01.2008
    def getApproximation( self, key, kind = 'Volume' ):
        iname, regionName, ig = key
##         print tregionName, aregionName, ig
#        print self.field.aps.apsPerGroup

        aps = self.field.aps
        geometries = aps.geometries
        ap = aps.apsPerGroup[ig]
        gKey = (iname, kind, regionName, ap.name)
        try:
            return ap, geometries[gKey]
        except KeyError:
            print gKey
            print geometries
            raise

    ##
    # c: 28.11.2006, r: 15.01.2008
    def getDataShapes( self, key, kind = 'Volume' ):
        iname, ig = key[0], key[-1]
        ap = self.field.aps.apsPerGroup[ig]
        if kind == 'Volume':
            shape = ap.getVDataShape( iname )
        else:
            regionName = key[1]
            shape = ap.getSDataShape( iname, regionName )

        return shape

    ##
    # 10.07.2007, c
    def getStateInRegion( self, region, igs = None, reshape = True ):
        nods = region.getFieldNodes( self.field, merge = True, igs = igs )
##         print nods, len( nods )
##         pause()
        eq = nm.empty( (len( nods ) * self.dpn,), dtype = nm.int32 )
        for idof in range( self.dpn ):
            eq[idof::self.dpn] = self.dpn * nods + idof + self.indx.start

        out = self.data[eq]
        if reshape:
            out.shape = (len( nods ), self.dpn)

        return out

    ##
    # 20.02.2007, c
    # 14.03.2007
    # 24.05.2007
    # 04.06.2007
    def extendData( self, data, nNod, val = None ):
        """Extend data (with value val) to cover whole domain."""
        cntVN = self.field.cntVN
        indx = cntVN[cntVN >= 0]

        if val is None:
            if data.shape[1] > 1: # Vector.
                val = nm.amin( nm.abs( data ) )
            else: # Scalar.
                val = nm.amin( data )

        extdata = nm.empty( (nNod, data.shape[1]), dtype = nm.float64 )
        extdata.fill( val )
        extdata[indx] = data[:indx.size]

        return extdata

## ##
## # 11.07.2006, c
## class FEVariable( Variable ):
##     """Finite element Variable
## field .. field description of variable (borrowed)
## """

# 24.09.2007, c
import sympy as s

##
# 25.09.2007, c
def createScalar( name, nEP ):
    vec = s.matrices.zeronm( nEP, 1 )
    for ip in range( nEP ):
        vec[ip,0] = '%s%d' % (name, ip)
    return vec

##
# 24.09.2007, c
def createVector( name, nEP, dim ):
    """ordering is DOF-by-DOF"""
    vec = s.matrices.zeronm( dim * nEP, 1 )
    for ii in range( dim ):
        for ip in range( nEP ):
            vec[nEP*ii+ip,0] = '%s%d%d' % (name, ii, ip)
    return vec

##
# 24.09.2007, c
def createScalarBase( name, nEP ):
    phi = s.matrices.zeronm( 1, nEP )
    for ip in range( nEP ):
        phi[0,ip] = '%s%d' % (name, ip)
    return phi

##
# 24.09.2007, c
# 25.09.2007
def createVectorBase( name, phic, dim ):
    nEP = phic.shape[1]
    phi = s.matrices.zeronm( dim, dim * nEP )
    indx = []
    for ii in range( dim ):
        phi[ii,nEP*ii:nEP*(ii+1)] = phic
        indx.append( ii )
    return phi, indx

##
# 24.09.2007, c
def createScalarBaseGrad( name, phic, dim ):
    nEP = phic.shape[1]
    gc = s.matrices.zeronm( dim, nEP )
    for ii in range( dim ):
        for ip in range( nEP ):
            gc[ii,ip] = '%s%d%d' % (name, ii, ip)
    return gc

##
# 24.09.2007, c
# 25.09.2007
def createVectorBaseGrad( name, gc, transpose = False ):
    dim, nEP = gc.shape
    g = s.matrices.zeronm( dim * dim, dim * nEP )
    indx = []
    if transpose:
        for ir in range( dim ):
            for ic in range( dim ):
                g[dim*ir+ic,nEP*ic:nEP*(ic+1)] = gc[ir,:]
                indx.append( (ic, ir) )
    else:
        for ir in range( dim ):
            for ic in range( dim ):
                g[dim*ir+ic,nEP*ir:nEP*(ir+1)] = gc[ic,:]
                indx.append( (ir, ic) )
    return g, indx

##
# 25.09.2007, c
def createUOperator( u, transpose = False ):
    dim = u.shape[0]
    opU = s.matrices.zeronm( dim * dim, dim )
    if transpose:
        for ir in range( dim ):
            for ic in range( dim ):
                opU[dim*ir+ic,ic] = u[ir]
    else:
        for ii in range( dim ):
            opU[dim*ii:dim*(ii+1),ii] = u
    return opU

##
# 24.09.2007, c
# 25.09.2007
def gradVectorToMatrix( name, gv ):
    dim2 = gv.shape[0]
    dim = int( s.sqrt( dim2 ) )
    gm = s.matrices.zeronm( dim, dim )
    for ir in range( dim ):
        for ic in range( dim ):
            gm[ir,ic] = gv[dim*ir+ic,0]
    return gm

##
# 24.09.2007, c
# 25.09.2007
def substituteContinuous( expr, names, u, phi ):
    pu = phi * u
    for ii in range( phi.lines ):
        expr = expr.subs( pu[ii,0], names[ii] )
    return expr

##
# 25.09.2007, c
def createVectorVarData( name, phi, vindx, g, gt, vgindx, u ):
    gu = g * u
    gum = gradVectorToMatrix( 'gum', gu )
    print 'g %s:\n' % name, gum

    gut = gt * u
    gutm = gradVectorToMatrix( 'gutm', gut )
    print 'gt %s:\n' % name, gutm

    pu = phi * u
    names = ['c%s%d' % (name, indx) for indx in vindx]
    cu = substituteContinuous( pu, names, u, phi )
    print 'continuous %s:\n' % name, cu

    gnames = ['cg%s%d_%d' % (name, indx[0], indx[1]) for indx in vgindx]
    cgu = substituteContinuous( gu, gnames, u, g )
    cgum = gradVectorToMatrix( 'gum', cgu )
    print 'continuous g %s:\n' % name, cgum

    cgut = substituteContinuous( gut, gnames, u, g )
    cgutm = gradVectorToMatrix( 'gutm', cgut )
    print 'continuous gt %s:\n' % name, cgutm

    opU = createUOperator( cu )
    print opU

    opUT = createUOperator( cu, transpose = True )
    print opUT

    out = {
        'g' : gu,
        'g_m' : gum,
        'q' : pu,
        'c' : cu,
        'cg' : cgu,
        'cg_m' : cgum,
        'cgt' : cgut,
        'cgt_m' : cgutm,
        'op' : opU,
        'opt' : opUT,
        'names' : names,
        'gnames' : gnames,
    }
    
    return out

##
# 25.09.2007, c
def createScalarVarData( name, phi, g, u ):
    gu = g * u

    pu = phi * u
    names = ['c%s' % name]
    cu = substituteContinuous( pu, names, u, phi )
    print 'continuous %s:\n' % name, cu

    gnames = ['cg%s_%d' % (name, ii) for ii in range( g.shape[0] )]
    cgu = substituteContinuous( gu, gnames, u, g )
    print 'continuous g %s:\n' % name, cgu

    opGU = createUOperator( cgu )
    print opGU

    out = {
        'g' : gu,
        'q' : pu,
        'c' : cu,
        'cg' : cgu,
        'gop' : opGU,
        'names' : names,
        'gnames' : gnames,
    }
    
    return out

##
# 25.09.2007, c
def main():
    nEP = 3
    dim = 2


    u = createVector( 'u', nEP, dim )
    v = createVector( 'v', nEP, dim )
    b = createVector( 'b', nEP, dim )
    p = createScalar( 'p', nEP )
    q = createScalar( 'q', nEP )
    r = createScalar( 'r', nEP )

    ## print u
    ## print v

    phic = createScalarBase( 'phic', nEP )
    phi, vindx = createVectorBase( 'phi', phic, dim )
    gc = createScalarBaseGrad( 'gc', phic, dim )
    g, vgindx = createVectorBaseGrad( 'g', gc )
    gt, aux = createVectorBaseGrad( 'gt', gc, transpose = True )

    ## print phi
    ## print phic
    ## print gc
    print g
    print gt

    ud = createVectorVarData( 'u', phi, vindx, g, gt, vgindx, u )
    vd = createVectorVarData( 'v', phi, vindx, g, gt, vgindx, v )
    bd = createVectorVarData( 'b', phi, vindx, g, gt, vgindx, b )
    pd = createScalarVarData( 'p', phic, gc, p )
    qd = createScalarVarData( 'q', phic, gc, q )
    rd = createScalarVarData( 'r', phic, gc, r )
    print ud.keys()

    assert bool( bd['op'].T * g == bd['opt'].T * gt )
    assert bool( bd['opt'].T * g == bd['op'].T * gt )
    assert bool( bd['cgt_m'] == bd['cg_m'].T )

    print '((b * grad) u), v)'
    form1 = vd['c'].T * bd['op'].T * ud['cg']
    form2 = vd['c'].T * bd['opt'].T * ud['cgt']
    print form1
    print form2
    print bool( form1 == form2 )

    print '((v * grad) u), b)'
    form1 = vd['c'].T * bd['op'].T * ud['cgt']
    form2 = vd['c'].T * bd['opt'].T * ud['cg']
    print form1
    print form2
    print bool( form1 == form2 )

    print '((u * grad) v), b)'
    form1 = vd['cgt'].T * bd['op'] * ud['c']
    form2 = vd['cg'].T * bd['opt'] * ud['c']
    print form1
    print form2
    print bool( form1 == form2 )

    print '((b * grad) v), u)'
    form1 = vd['cg'].T * bd['op'] * ud['c']
    form2 = vd['cgt'].T * bd['opt'] * ud['c']
    print form1
    print form2
    print bool( form1 == form2 )

    print '((v * grad) b), u)'
    form1 = vd['c'].T * bd['cgt_m'] * ud['c']
    form2 = vd['c'].T * bd['cg_m'].T * ud['c']
    print form1
    print form2
    print bool( form1 == form2 )

    print '((b * grad) u), (b * grad) v)'
    form1 = vd['cg'].T * bd['op'] * bd['op'].T * ud['cg']
    print form1

    print '((u * grad) b), (b * grad) v)'
    form1 = vd['cg'].T * bd['op'] * bd['cg_m'] * ud['c']
    print form1

    print '(grad p, (b * grad) v)'
    form1 = vd['cg'].T * bd['op'] * pd['cg']
    print form1
    
    print '(grad q, (b * grad) u)'
    form1 = qd['cg'].T * bd['op'].T * ud['cg']
    print form1

    print '(grad q, (u * grad) b)'
    form1 = qd['cg'].T * bd['cg_m'] * ud['c']
    print form1

    print '(grad r, (u * grad) v)'
    form1 = vd['cgt'].T * rd['gop'] * ud['c']
    print form1
    
    return ud, vd, bd, pd, qd, rd

if __name__ == '__main__':
    ud, vd, bd, pd, qd, rd = main()

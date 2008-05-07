from base import *
import inspect

##
# 30.05.2007, c
class TestCommon( Struct ):

    ##
    # 16.07.2007, c
    def getNumber( self ):
        methods = inspect.getmembers( self, inspect.ismethod )
        tests = [ii for ii in methods
                 if (len( ii[0] ) > 5) and ii[0][:5] == 'test_']
        return len( tests )

    ##
    # c: 30.05.2007, r: 05.02.2008
    def run( self, debug = False ):
        ok = True
        nFail = 0

        methods = inspect.getmembers( self, inspect.ismethod )
        if hasattr( self, 'tests' ):
            dmethods = {}
            for key, method in methods:
                dmethods[key] = method
            tests = [(ii, dmethods[ii]) for ii in self.tests]
            print tests
        else:
            tests = [ii for ii in methods
                     if (len( ii[0] ) > 5) and ii[0][:5] == 'test_']

            
        for testName, testMethod in tests:
            aux = ' %s: ' % testName

            try:
                ret = testMethod()
            except:
                if debug:
                    raise
                ret = False
                
            if not ret:
                aux = '---' + aux + 'failed!'
                nFail += 1
                ok = False
            else:
                aux = '+++' + aux + 'ok'

            print aux

        return ok, nFail, len( tests )

    ##
    # c: 31.05.2007, r: 02.05.2008
    def report( *argc ):
        """All tests should print via this function."""
        format = '...' + ' %s' * len( argc )
        msg =  format % argc
        print msg
    report = staticmethod( report )

    ##
    # 30.05.2007, c
    def evalCoorExpression( expression, coor ):

        x = coor[:,0]
        y = coor[:,1]
        if coor.shape[1] == 3:
            z = coor[:,2]
        else:
            z = None

        coorDict = {'x' : x, 'y' : y, 'z' : z}
        out = eval( expression, globals(), coorDict )
        
        return out
    evalCoorExpression = staticmethod( evalCoorExpression )

    ##
    # c: 30.05.2007, r: 07.05.2008
    def compareVectors( vec1, vec2, allowedError = 1e-8,
                        label1 = 'vec1', label2 = 'vec2', norm = None ):

        diffNorm = nla.norm( vec1 - vec2, ord = norm )
        TestCommon.report( '||%s - %s||: %e' % (label1, label2, diffNorm) )
        if diffNorm > allowedError:
            return False
        else:
            return True
    compareVectors = staticmethod( compareVectors )

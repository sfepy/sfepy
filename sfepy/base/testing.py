from base import *
import inspect

##
# 30.05.2007, c
class TestCommon( Struct ):

    ##
    # 16.07.2007, c
    def get_number( self ):
        methods = inspect.getmembers( self, inspect.ismethod )
        tests = [ii for ii in methods
                 if (len( ii[0] ) > 5) and ii[0][:5] == 'test_']
        return len( tests )

    ##
    # c: 30.05.2007, r: 05.02.2008
    def run( self, debug = False ):
        ok = True
        n_fail = 0

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

            
        for test_name, test_method in tests:
            aux = ' %s: ' % test_name

            try:
                ret = test_method()
            except:
                if debug:
                    raise
                ret = False
                
            if not ret:
                aux = '---' + aux + 'failed!'
                n_fail += 1
                ok = False
            else:
                aux = '+++' + aux + 'ok'

            print aux

        return ok, n_fail, len( tests )

    ##
    # c: 31.05.2007, r: 02.05.2008
    def report( *argc ):
        """All tests should print via this function."""
        format = '...' + ' %s' * len( argc )
        msg =  format % argc
        print msg
    report = staticmethod( report )

    ##
    # c: 30.05.2007, r: 09.05.2008
    def eval_coor_expression( expression, coor ):

        x = coor[:,0]
        y = coor[:,1]
        if coor.shape[1] == 3:
            z = coor[:,2]
        else:
            z = None

        env = {'x' : x, 'y' : y, 'z' : z}
        out = eval( expression, nm.__dict__, env )
        
        return out
    eval_coor_expression = staticmethod( eval_coor_expression )

    ##
    # c: 30.05.2007, r: 07.05.2008
    def compare_vectors( vec1, vec2, allowed_error = 1e-8,
                        label1 = 'vec1', label2 = 'vec2', norm = None ):

        diff_norm = nla.norm( vec1 - vec2, ord = norm )
        TestCommon.report( '||%s - %s||: %e' % (label1, label2, diff_norm) )
        if diff_norm > allowed_error:
            return False
        else:
            return True
    compare_vectors = staticmethod( compare_vectors )

from sfepy.base.base import assert_
from sfepy.base.testing import TestCommon

##
# 28.08.2007, c
class Test( TestCommon ):

    ##
    # 28.08.2007, c
    def from_conf( conf, options ):
        return Test( conf = conf, options = options )
    from_conf = staticmethod( from_conf )

    ##
    # 28.08.2007, c
    def test_struct_add( self ):
        from sfepy.base.base import Struct
        from copy import deepcopy

        a = Struct( f1 = 0,
                    f2 = [1, 2, 3],
                    f3 = Struct( ff = 'abc' ),
                    f4 = 3.14 )
        a0 = deepcopy( a )
        b = Struct( f1 = 5,
                    f2 = [1],
                    f3 = Struct( ff = '', gg = 123 ),
                    f5 = 'new one' )
        c = a + b

        assert_( c.f1 == 0 )
        assert_( c.f2 == [1, 2, 3] )
        assert_( c.f3.ff == 'abc' )
        assert_( c.f3.gg == 123 )
        assert_( c.f4 == 3.14 )
        assert_( c.f5 == 'new one' )

        assert_( a.f1 == a0.f1 )
        assert_( a.f2 == a0.f2 )
        assert_( a.f3.ff == a0.f3.ff )
        assert_( a.f4 == a0.f4 )

        return True

    ##
    # 28.08.2007, c
    def test_struct_i_add( self ):
        from sfepy.base.base import Struct

        a = Struct( f1 = 0,
                    f2 = [1, 2, 3],
                    f3 = Struct( ff = 'abc' ) )
        b = Struct( f1 = 5,
                    f2 = [1],
                    f3 = Struct( ff = '', gg = 123 ),
                    f4 = 'new one' )
        a += b

        assert_( a.f1 == 0 )
        assert_( a.f2 == [1, 2, 3] )
        assert_( a.f3.ff == 'abc' )
        assert_( a.f3.gg == 123 )
        assert_( a.f4 == 'new one' )

        return True

    def test_verbose_output(self):
        import StringIO
        from sfepy.base.base import Output, goptions

        fd = StringIO.StringIO()

        output = Output('test', filename=fd)

        output('test1')
        goptions['verbose'] = False
        output('test2')
        goptions['verbose'] = 1
        output('test3')

        _ok1 = goptions['verbose'] == True
        _ok2 = fd.getvalue() == 'test test1\ntest test3\n'

        fd.close()

        ok = _ok1 and _ok2

        return ok

    def test_resolve_deps(self):
        from sfepy.base.resolve_deps import resolve

        deps = {
            'a' : ['a', 'b'],
            'b' : ['a', 'b'],
            'c' : ['b', 'c', 'd', 'e'],
            'd' : ['c', 'e'],
            'e' : ['d', 'e'],
            'f' : ['e', 'f', 'g'],
            'g' : ['g'],
        }
        order = resolve(deps)
        ok1 = order == [['g'], ['a', 'b'], ['c', 'd', 'e'], ['f']]

        deps = {
            'a' : ['b'],
            'b' : ['c'],
            'c' : ['a'],
        }
        order = resolve(deps)
        ok2 = order == [['a', 'b', 'c']]

        return ok1 and ok2

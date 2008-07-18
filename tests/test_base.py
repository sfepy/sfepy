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

        assert c.f1 == 0
        assert c.f2 == [1, 2, 3]
        assert c.f3.ff == 'abc'
        assert c.f3.gg == 123
        assert c.f4 == 3.14
        assert c.f5 == 'new one'

        assert a.f1 == a0.f1
        assert a.f2 == a0.f2
        assert a.f3.ff == a0.f3.ff
        assert a.f4 == a0.f4

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

        assert a.f1 == 0
        assert a.f2 == [1, 2, 3]
        assert a.f3.ff == 'abc'
        assert a.f3.gg == 123
        assert a.f4 == 'new one'

        return True

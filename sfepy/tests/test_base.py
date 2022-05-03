from sfepy.base.base import assert_

def test_struct_add():
    from sfepy.base.base import Struct
    from copy import deepcopy

    a = Struct(f1=0,
               f2=[1, 2, 3],
               f3=Struct(ff='abc'),
               f4=3.14)
    a0 = deepcopy(a)
    b = Struct(f1=5,
               f2=[1],
               f3=Struct(ff='', gg=123),
               f5='new one')
    c = a + b

    assert_(c.f1 == 0)
    assert_(c.f2 == [1, 2, 3])
    assert_(c.f3.ff == 'abc')
    assert_(c.f3.gg == 123)
    assert_(c.f4 == 3.14)
    assert_(c.f5 == 'new one')

    assert_(a.f1 == a0.f1)
    assert_(a.f2 == a0.f2)
    assert_(a.f3.ff == a0.f3.ff)
    assert_(a.f4 == a0.f4)

def test_struct_i_add():
    from sfepy.base.base import Struct

    a = Struct(f1=0,
               f2=[1, 2, 3],
               f3=Struct(ff='abc'))
    b = Struct(f1=5,
               f2=[1],
               f3=Struct(ff='', gg=123),
               f4='new one')
    a += b

    assert_(a.f1 == 0)
    assert_(a.f2 == [1, 2, 3])
    assert_(a.f3.ff == 'abc')
    assert_(a.f3.gg == 123)
    assert_(a.f4 == 'new one')

def test_container_add():
    from sfepy.base.base import Struct, Container

    a = Struct(name='a')
    b = Struct(name='b')

    c1 = Container()

    c1 = c1 + c1
    assert_(c1.names == [])

    c1 += Container([a, b])
    assert_(c1.names == ['a', 'b'])

    c2 = c1 + c1
    assert_(c2.names == 2 * ['a', 'b'])

    c2 += c2
    assert_(c2.names == 4 * ['a', 'b'])

def test_verbose_output():
    from io import StringIO

    from sfepy.base.base import Output, goptions

    fd = StringIO()

    output = Output('test', filename=fd)

    output('test1')
    goptions['verbose'] = False
    output('test2')
    goptions['verbose'] = 1
    output('test3')

    _ok1 = bool(goptions['verbose'])
    _ok2 = fd.getvalue() == 'test test1\ntest test3\n'

    fd.close()

    ok = _ok1 and _ok2

    assert_(ok)

def test_resolve_deps():
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

    assert_(ok1 and ok2)

def test_parse_conf():
    from sfepy.base.parse_conf import list_dict
    from sfepy.base.conf import dict_from_string as parse_dict
    ld = list_dict()

    def parse(x):
        out = ld.parseString(x, True)[0]
        return out

    assert_(parse("1,2") == ([1,2],{}))
    assert_(parse("1,[2]") == ([1,[2]],{}))
    assert_(parse("1,[2,4],c=3") == ([1,[2,4]],{'c':3}))
    assert_(parse("1,(2),c:3,uu=7") == ([1,(2,)],{'c':3,'uu':7}))
    assert_(parse("'long string ([\"',(2,5),c:3") ==
                 (['long string (["',(2,5)],{'c':3}))


    assert_(parse_dict('') == {})
    assert_(parse_dict('a:[]') == {'a':[]})
    assert_(parse_dict('a:{}') == {'a':{}})
    assert_(parse_dict('1:2,a:{},3:4') == {1:2,'a':{},3:4})

    assert_(parse('') == ([],{}))
    assert_(parse('[ ]') == ([[]],{}))
    assert_(parse('[]') == ([[]],{}))
    assert_(parse('[[]]') == ([[[]]],{}))
    assert_(parse('[[[]]]')==([[[[]]]],{}))
    assert_(parse('a,{},[],None,True,False,"False"') ==
                     (['a',{},[],None,True,False,"False"],{}))

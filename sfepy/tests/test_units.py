from sfepy.base.base import assert_
import sfepy.base.testing as tst

def _cmp(s1, s2):
    s1 = s1.split()
    s2 = s2.split()
    v1, t1 = float(s1[0]), s1[1:]
    v2, t2 = float(s2[0]), s2[1:]
    return (abs(v1 - v2) < (1e-15 * abs(v1))) and (t1 == t2)

def test_units():
    from sfepy.mechanics.units import Unit, Quantity, sm

    if sm is None:
        tst.report('cannot import sympy, skipping')
        return

    units = ['m', 's', 'kg', 'C']
    tst.report('units:', units)
    unit_set = [Unit(key) for key in units]

    q1 = Quantity('stress', unit_set)
    tst.report(q1.name, ':', q1())
    assert_(_cmp(q1(), '1.0 Pa'))
    assert_(_cmp(q1('c'), '100.0 cPa'))

    q2 = Quantity('force', unit_set)
    tst.report(q2.name, ':', q2())
    assert_(_cmp(q2(), '1.0 Newton'))
    assert_(_cmp(q2('d'), '0.1 dNewton'))

    q3 = Quantity('energy', unit_set)
    tst.report(q3.name, ':', q3())
    assert_(_cmp(q3(), '1.0 J'))
    assert_(_cmp(q3('mu'), '1000000.0 muJ'))

    units = ['mm', 's', 'g', 'C']
    tst.report('units:', units)
    unit_set = [Unit(key) for key in units]

    q1 = Quantity('stress', unit_set)
    tst.report(q1.name, ':', q1())
    assert_(_cmp(q1(), '1.0 Pa'))

    q2 = Quantity('force', unit_set)
    tst.report(q2.name, ':', q2())
    assert_(_cmp(q2(), '1.0 muNewton'))

    q3 = Quantity('energy', unit_set)
    tst.report(q3.name, ':', q3())
    assert_(_cmp(q3(), '1.0 nJ'))

    units = ['cm', 'ms', 'kg', 'kC']
    tst.report('units:', units)
    unit_set = [Unit(key) for key in units]

    q1 = Quantity('stress', unit_set)
    tst.report(q1.name, ':', q1())
    assert_(_cmp(q1(), '0.1 GPa'))

    q2 = Quantity('force', unit_set)
    tst.report(q2.name, ':', q2())
    assert_(_cmp(q2(), '10.0 kNewton'))

    q3 = Quantity('energy', unit_set)
    tst.report(q3.name, ':', q3())
    assert_(_cmp(q3(), '0.1 kJ'))

    q4 = Quantity('thermal_expandability', unit_set)
    tst.report(q4.name, ':', q4())
    assert_(_cmp(q4(), '0.1 MPa / C'))

    assert_(_cmp(q4('G'), '0.0001 GPa / C'))
    assert_(_cmp(q4('M'), '0.1 MPa / C'))
    assert_(_cmp(q4('k'), '100.0 kPa / C'))
    assert_(_cmp(q4('d'), '10000.0 dPa / C'))
    assert_(_cmp(q4(''), '100000.0 Pa / C'))

    units = ['m', 's', 'g', 'C']
    tst.report('units:', units)
    unit_set = [Unit(key) for key in units]

    q4 = Quantity('thermal_expandability', unit_set)
    tst.report(q4.name, ':', q4())
    assert_(_cmp(q4(), '1.0 mPa / C'))

    assert_(_cmp(q4('k'), str(0.000001) + ' kPa / C'))
    assert_(_cmp(q4('d'), '0.0001 dPa / C'))
    assert_(_cmp(q4(''), '0.001 Pa / C'))
    assert_(_cmp(q4('c'), '0.1 cPa / C'))
    assert_(_cmp(q4('m'), '1.0 mPa / C'))
    assert_(_cmp(q4('mu'), '1000.0 muPa / C'))
    assert_(_cmp(q4('n'), '1000000.0 nPa / C'))

def test_consistent_sets():
    from sfepy.mechanics.units import get_consistent_unit_set, sm

    if sm is None:
        tst.report('cannot import sympy, skipping')
        return True

    u_sets = {
        ('m', 's', 'kg', 'C') : {'force' : '1.0 Newton',
                                 'stress' : '1.0 Pa',
                                 'energy' : '1.0 J',
                                 'thermal_expandability' : '1.0 Pa / C'},
        ('mm', 's', 'kg', 'C') : {'force' : '1.0 mNewton',
                                  'stress' : '1.0 kPa',
                                  'energy' : '1.0 muJ',
                                 'thermal_expandability' : '1.0 kPa / C'},
        ('mm', 's', 'g', 'C') : {'force' : '1.0 muNewton',
                                 'stress' : '1.0 Pa',
                                 'energy' : '1.0 nJ',
                                 'thermal_expandability' : '1.0 Pa / C'},
    }

    ok = True
    for unit_set, true_derived_units in u_sets.items():
        tst.report('units:', unit_set)
        derived_units = get_consistent_unit_set(*unit_set)

        for key, true_val in true_derived_units.items():
            val = derived_units[key]
            _ok = _cmp(true_val, val)
            tst.report('%s: %s == %s -> %s' % (key, true_val, val, _ok))

            ok = ok and _ok

    assert_(ok)

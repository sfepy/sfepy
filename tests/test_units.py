from sfepy.base.testing import TestCommon, assert_, debug

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        return Test(conf=conf, options=options)

    def test_units(self):
        from sfepy.mechanics.units import Unit, Quantity

        units = ['m', 's', 'kg', 'C']
        self.report('units:', units)
        unit_set = [Unit(key) for key in units]
 
        q1 = Quantity('stress', unit_set)
        self.report(q1.name, ':', q1())
        assert_(q1() == '1.0 Pa')
        assert_(q1('c') == '100.0 cPa')

        q2 = Quantity('force', unit_set)
        self.report(q2.name, ':', q2())
        assert_(q2() == '1.0 Newton')
        assert_(q2('d') == '0.1 dNewton')

        q3 = Quantity('energy', unit_set)
        self.report(q3.name, ':', q3())
        assert_(q3() == '1.0 J')
        assert_(q3('mu') == '1000000.0 muJ')

        units = ['mm', 's', 'g', 'C']
        self.report('units:', units)
        unit_set = [Unit(key) for key in units]
 
        q1 = Quantity('stress', unit_set)
        self.report(q1.name, ':', q1())
        assert_(q1() == '1.0 Pa')

        q2 = Quantity('force', unit_set)
        self.report(q2.name, ':', q2())
        assert_(q2() == '1.0 muNewton')

        q3 = Quantity('energy', unit_set)
        self.report(q3.name, ':', q3())
        assert_(q3() == '1.0 nJ')

        units = ['cm', 'ms', 'kg', 'kC']
        self.report('units:', units)
        unit_set = [Unit(key) for key in units]
 
        q1 = Quantity('stress', unit_set)
        self.report(q1.name, ':', q1())
        assert_(q1() == '0.1 GPa')

        q2 = Quantity('force', unit_set)
        self.report(q2.name, ':', q2())
        assert_(q2() == '10.0 kNewton')

        q3 = Quantity('energy', unit_set)
        self.report(q3.name, ':', q3())
        assert_(q3() == '0.1 kJ')

        q4 = Quantity('thermal_expandability', unit_set)
        self.report(q4.name, ':', q4())
        assert_(q4() == '0.1 MPa / C')

        assert_(q4('G') == '0.0001 GPa / C')
        assert_(q4('M') == '0.1 MPa / C')
        assert_(q4('k') == '100.0 kPa / C')
        assert_(q4('d') == '10000.0 dPa / C')
        assert_(q4('') == '100000.0 Pa / C')

        units = ['m', 's', 'g', 'C']
        self.report('units:', units)
        unit_set = [Unit(key) for key in units]

        q4 = Quantity('thermal_expandability', unit_set)
        self.report(q4.name, ':', q4())
        assert_(q4() == '1.0 mPa / C')

        assert_(q4('k') == str(0.000001) + ' kPa / C')
        assert_(q4('d') == '0.0001 dPa / C')
        assert_(q4('') == '0.001 Pa / C')
        assert_(q4('c') == '0.1 cPa / C')
        assert_(q4('m') == '1.0 mPa / C')
        assert_(q4('mu') == '1000.0 muPa / C')
        assert_(q4('n') == '1000000.0 nPa / C')

        return True


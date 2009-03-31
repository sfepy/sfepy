import re

from base import Struct, IndexedStruct, dict_to_struct, pause, output, copy,\
     import_file, assert_, get_default
from reader import Reader

_required = ['filename_mesh', 'field_[0-9]+|fields',
             'ebc_[0-9]+|ebcs', 'fe', 'equations',
             'region_[0-9]+|regions', 'variable_[0-9]+|variables',
             'material_[0-9]+|materials', 'integral_[0-9]+|integrals',
             'solver_[0-9]+|solvers']
_other = ['epbc_[0-9]+|epbcs', 'lcbc_[0-9]+|lcbcs', 'nbc_[0-9]+|nbcs',
          'ic_[0-9]+|ics', 'options']

##
# c: 19.02.2008, r: 19.02.2008
def get_standard_keywords():
    return copy( _required ), copy( _other )

##
# c: 10.04.2008, r: 10.04.2008
def tuple_to_conf( name, vals, order ):
    conf = Struct( name = name )
    for ii, key in enumerate( order ):
        setattr( conf, key, vals[ii] )
    return conf

##
# Short syntax: key is suffixed with '__<number>' to prevent collisions with
# long syntax keys -> both cases can be used in a single input.
def transform_variables( adict ):
    d2 = {}
    for ii, (key, conf) in enumerate( adict.iteritems() ):
        if isinstance( conf, tuple ):
            c2 = tuple_to_conf( key, conf, ['kind', 'field'] )
            if len( conf ) >= 3:
                kind = c2.kind.split()[0]
                if kind == 'unknown':
                    c2.order = conf[2]
                elif kind == 'test':
                    c2.dual = conf[2]
                elif kind == 'parameter':
                    c2.like = conf[2]
                if len( conf ) == 4:
                    c2.history = conf[3]
            d2['variable_%s__%d' % (c2.name, ii)] = c2
        else:
            c2 = transform_to_struct_1( conf )
            d2['variable_'+c2.name] = c2
    return d2

def transform_conditions( adict, prefix ):
    d2 = {}
    for ii, (key, conf) in enumerate( adict.iteritems() ):
        if isinstance( conf, tuple ):
            c2 = tuple_to_conf( key, conf, ['region', 'dofs'] )
            if len( conf ) == 3:
                c2.filename = conf[2]
            d2['%s_%s__%d' % (prefix, c2.name, ii)] = c2
        else:
            c2 = transform_to_struct_1( conf )
            d2['%s_%s' % (prefix, c2.name)] = c2
    return d2
    
def transform_ebcs( adict ):
    return transform_conditions( adict, 'ebc' )

def transform_ics( adict ):
    return transform_conditions( adict, 'ic' )

def transform_lcbcs( adict ):
    return transform_conditions( adict, 'lcbc' )

def transform_regions( adict ):
    d2 = {}
    for ii, (key, conf) in enumerate( adict.iteritems() ):
        if isinstance( conf, tuple ):
            c2 = tuple_to_conf( key, conf, ['select', 'flags'] )
            for flag, val in c2.flags.iteritems():
                setattr( c2, flag, val )
            delattr( c2, 'flags' )
            d2['region_%s__%d' % (c2.name, ii)] = c2
        else:
            c2 = transform_to_struct_1( conf )
            d2['region_'+c2.name] = c2
    return d2

def transform_integrals( adict ):
    d2 = {}
    for ii, (key, conf) in enumerate( adict.iteritems() ):
        if isinstance( conf, tuple ):
            c2 = tuple_to_conf( key, conf, ['kind', 'quadrature'] )
            d2['integral_%s__%d' % (c2.name, ii)] = c2
        else:
            c2 = transform_to_struct_1( conf )
            d2['integral_'+c2.name] = c2
    return d2

def transform_fields( adict ):
    d2 = {}
    for ii, (key, conf) in enumerate( adict.iteritems() ):
        if isinstance( conf, tuple ):
            c2 = tuple_to_conf( key, conf, ['dim', 'domain','bases'] )
            d2['field_%s__%d' % (c2.name, ii)] = c2
        else:
            c2 = transform_to_struct_1( conf )
            d2['field_'+c2.name] = c2
    return d2

def transform_materials( adict ):
    d2 = {}
    for ii, (key, conf) in enumerate( adict.iteritems() ):
        if isinstance( conf, tuple ):
            c2 = tuple_to_conf( key, conf, ['mode', 'region', 'params'] )
            dpar = {}
            dpar['name'] = key
            dpar['mode'] = c2.mode
            dpar['region'] = c2.region
            dpar.update( c2.params )
            d2['material_%s__%d' % (c2.name, ii)] = dpar
        else:
            d2['material_'+conf['name']] = conf
    return d2

def transform_solvers( adict ):
    d2 = {}
    for ii, (key, conf) in enumerate( adict.iteritems() ):
        if isinstance( conf, tuple ):
            c2 = tuple_to_conf( key, conf, ['kind','params'] )
            for param, val in c2.params.iteritems():
                setattr( c2, param, val )
            delattr( c2, 'params' )
            d2['solvers_%s__%d' % (c2.name, ii)] = c2
        else:
            c2 = transform_to_struct_1( conf )
            d2['solvers_'+c2.name] = c2
    return d2

##
# c: 20.06.2007, r: 18.02.2008
def transform_to_struct_1( adict ):
    return dict_to_struct( adict, flag = (1,) )
def transform_to_i_struct_1( adict ):
    return dict_to_struct( adict, flag = (1,), constructor = IndexedStruct )
def transform_to_struct_01( adict ):
    return dict_to_struct( adict, flag = (0,1) )
def transform_to_struct_10( adict ):
    return dict_to_struct( adict, flag = (1,0) )

transforms = {
    'options'   : transform_to_i_struct_1,
    'solvers'   : transform_solvers,
    'integrals' : transform_integrals,
    'opt'       : transform_to_struct_1,
    'fe'        : transform_to_struct_1,
    'regions'   : transform_regions,
    'shape_opt' : transform_to_struct_10,
    'fields'    : transform_fields,
    'variables' : transform_variables,
    'ebcs'      : transform_ebcs,
    'epbcs'     : transform_to_struct_01,
    'nbcs'      : transform_to_struct_01,
    'lcbcs'     : transform_lcbcs,
    'ics'       : transform_ics,
    'materials' : transform_materials,
}

##
# 27.10.2005, c
class ProblemConf( Struct ):
    """
    Problem configuration, corresponding to an input (problem description
    file). It validates the input using lists of required and other keywords
    that have to/can appear in the input. Default keyword lists can be obtained
    by sfepy.base.conf.get_standard_keywords().

    ProblemConf instance is used to construct a ProblemDefinition instance via
    ProblemDefinition.from_conf( conf ).
    """
    ##
    # c: 25.07.2006, r: 10.07.2008
    def from_file( filename, required = None, other = None ):
        """
        Loads the problem definition from a file.

        The filename can either contain plain definitions, or it can contain
        the define() function, in which case it will be called to return the
        input definitions.

        The job of the define() function is to return a dictionary of
        parameters. How the dictionary is constructed is not our business, but
        the usual way is to simply have a function define() along these lines
        in the input file:

            def define():
                options = {
                    'save_eig_vectors' : None,
                    'eigen_solver' : 'eigen1',
                }
                region_2 = {
                    'name' : 'Surface',
                    'select' : 'nodes of surface',
                }
                return locals()

        """
        funmod = import_file( filename )
        obj = ProblemConf()
        if "define" in funmod.__dict__:
            define_dict = funmod.__dict__["define"]()
        else:
            define_dict = funmod.__dict__

        obj.__dict__.update( define_dict )

        obj.setup( define_dict, funmod, filename, required, other )

        return obj
    from_file = staticmethod( from_file )

    def from_module( module, required = None, other = None ):
        obj = ProblemConf()
        obj.__dict__.update( module.__dict__ )

        obj.setup( funmod = module, required = required, other = other )

        return obj
    from_module = staticmethod( from_module )

    def from_dict( dict_, funmod, required = None, other = None ):
        obj = ProblemConf()
        obj.__dict__.update( dict_ )

        obj.setup( funmod = funmod, required = required, other = other )

        return obj
    from_dict = staticmethod( from_dict )

    def setup( self, define_dict = None, funmod = None, filename = None,
               required = None, other = None ):

        define_dict = get_default( define_dict, self.__dict__ )

        self._filename = filename

        other_missing = self.validate( required = required, other = other )

        self.transform_input_trivial()
        self._raw = {}
        for key, val in define_dict.iteritems():
            if isinstance( val, dict ):
                self._raw[key] = copy( val )

        self.transform_input()
        self.funmod = funmod

    ##
    # 27.10.2005, c
    # 19.09.2006
    # 05.06.2007
    def _validate_helper( self, items, but_nots ):
        keys = self.__dict__.keys()
        left_over = keys[:]
        if but_nots is not None:
            for item in but_nots:
                match = re.compile( '^' + item + '$' ).match
                for key in keys:
                    if match( key ):
                        left_over.remove( key )

        missing = []
        if items is not None:
            for item in items:
                found = False
                match = re.compile( '^' + item + '$' ).match
                for key in keys:
                    if match( key ):
                        found = True
                        left_over.remove( key )
                if not found:
                    missing.append( item )
        return left_over, missing

    ##
    # c: 27.10.2005, r: 11.07.2008
    def validate( self, required = None, other = None ):
        required_left_over, required_missing \
                          = self._validate_helper( required, other )
        other_left_over, other_missing \
                       = self._validate_helper( other, required )

        assert_( required_left_over == other_left_over )

        err = False
        if required_missing:
            err = True
            output( 'error: required missing:', required_missing )

        if other_left_over:
            output( 'left over:', other_left_over )

        if err:
            raise ValueError

        return other_missing

    ##
    # c: 31.10.2005, r: 10.07.2008
    def transform_input_trivial( self ):
        """Trivial input transformations."""

        ##
        # Unordered inputs.
        tr_list = ['([a-zA-Z0-9]+)_[0-9]+']
        # Keywords not in 'required', but needed even empty (e.g. for run_tests).
        for key in transforms.keys():
            if not self.__dict__.has_key( key ):
                self.__dict__[key] = {}

        keys = self.__dict__.keys()
        for item in tr_list:
            match = re.compile( item ).match
            for key in keys:
                obj = match( key )
                if obj:
                    new = obj.group( 1 ) + 's'
                    result = {key : self.__dict__[key]}
                    try:
                        self.__dict__[new].update( result )
                    except:
                        self.__dict__[new] = result
                        
                    del self.__dict__[key]

    def transform_input( self ):
        keys = self.__dict__.keys()
        for key, transform in transforms.iteritems():
            if not key in keys: continue
            self.__dict__[key] = transform( self.__dict__[key] )

    def get_raw( self, key = None ):
        if key is None:
            return self._raw
        else:
            return self._raw[key]

    def edit( self, key, newval ):
        self.__dict__[key] = transforms[key]( newval )

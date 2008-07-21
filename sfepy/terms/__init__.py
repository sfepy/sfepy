import extmods
import terms
from termsBasic import *
from termsLinElasticity import *
from termsLaplace import *
from termsPoint import *
from termsSurface import *
from termsVolume import *
from termsMass import *
from termsNavierStokes import *
from termsBiot import *

try:
    from termsAdjointNavierStokes import *
except:
    pass

try:
    from termsHDPM import *
    from cachesHDPM import *
except:
    pass

try:
    from cachesBasic import *
except:
    from sfepy.base.base import output
    output( 'warning: sfepy extension modules are not compiled!' )
    output( 'type "make"' )

##
# 15.11.2005, c
# 16.11.2005
# 28.11.2005
# 24.11.2006
# 27.02.2007
var_dict = vars().items()
term_table = {}
cache_table = {}
is_term = re.compile( '[a-zA-Z_0-9]+Term$' ).match
is_cache = re.compile( '[a-zA-Z_0-9]+DataCache$' ).match

for key, var in var_dict:
#    print key, var
    if is_term( key ):
        term_table[var.name] = var
    elif is_cache( key ):
        cache_table[var.name] = var
del var_dict
## print term_table.keys()
## print cache_table.keys()
## pause()

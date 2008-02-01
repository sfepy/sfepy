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

try:
    from termsAdjointNavierStokes import *
except:
    pass

try:
    from termsHDPM import *
    from cachesHDPM import *
except:
    pass

from cachesBasic import *

##
# 15.11.2005, c
# 16.11.2005
# 28.11.2005
# 24.11.2006
# 27.02.2007
varDict = vars().items()
termTable = {}
cacheTable = {}
isTerm = re.compile( '[a-zA-Z_0-9]+Term$' ).match
isCache = re.compile( '[a-zA-Z_0-9]+DataCache$' ).match

for key, var in varDict:
#    print key, var
    if isTerm( key ):
        termTable[var.name] = var
    elif isCache( key ):
        cacheTable[var.name] = var
del varDict
## print termTable.keys()
## print cacheTable.keys()
## pause()

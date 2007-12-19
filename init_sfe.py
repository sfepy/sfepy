##
# Path and module setup.
# 18.01.2005
# 13.10.2005
import sys
import os.path as op

install_dir = op.dirname( __file__ )

#sys.path.append( install_dir + "/sfelib" )
#sys.path.append( install_dir + "/eldesc" )

##
# Taken from pyunit.
## import __builtin__
## class RollbackImporter:
##     def __init__(self):
##         "Creates an instance and installs as the global importer"
##         self.previousModules = sys.modules.copy()
##         self.realImport = __builtin__.__import__
##         __builtin__.__import__ = self._import
##         self.newModules = {}
        
##     def _import(self, name, globals=None, locals=None, fromlist=[]):
##         result = apply(self.realImport, (name, globals, locals, fromlist))
##         self.newModules[name] = 1
##         return result
        
##     def uninstall(self):
##         for modname in self.newModules.keys():
##             if not self.previousModules.has_key(modname):
##                 # Force reload when modname next imported
##                 del(sys.modules[modname])
##         __builtin__.__import__ = self.realImport

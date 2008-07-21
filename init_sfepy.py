##
# Path and module setup.
# 18.01.2005
# 13.10.2005
import sys
import os.path as op

install_dir = op.dirname( __file__ )

##
# Taken from pyunit.
## import __builtin__
## class RollbackImporter:
##     def __init__(self):
##         "Creates an instance and installs as the global importer"
##         self.previous_modules = sys.modules.copy()
##         self.real_import = __builtin__.__import__
##         __builtin__.__import__ = self._import
##         self.new_modules = {}
        
##     def _import(self, name, globals=None, locals=None, fromlist=[]):
##         result = apply(self.real_import, (name, globals, locals, fromlist))
##         self.new_modules[name] = 1
##         return result
        
##     def uninstall(self):
##         for modname in self.new_modules.keys():
##             if not self.previous_modules.has_key(modname):
##                 # Force reload when modname next imported
##                 del(sys.modules[modname])
##         __builtin__.__import__ = self.real_import

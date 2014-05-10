import os
from cffi import FFI
from numpy import *

try:
   from pydb import debugger

   ## Also add an exception hook.
   import pydb, sys
   sys.excepthook = pydb.exception_hook

except ImportError:
   import pdb
   def debugger():
       pdb.set_trace()
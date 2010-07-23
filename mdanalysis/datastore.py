from tables import *
import numpy
import os
import inspect

_tables = []

#
# PyTables table definitions
#

class ExampleTable(IsDescription):
     name      = StringCol(16)   # 16-character String
     signed_64bit  = Int64Col()      # Signed 64-bit integer
     short_int  = UInt16Col()     # Unsigned short integer
     unsigned_byte  = UInt8Col()      # unsigned byte
     integer    = Int32Col()      # 32-bit integer
     float_sp  = Float32Col()    # float  (single-precision)
     double_dp    = Float64Col()    # double (double-precision)

#
# Initialize the table file
#

def create(filename, module=self, title="datastore"):
    """ 
    Create the pytables data file at the given filename. It will not run if the file already exists.
    For each class defined in this file which inherits IsDescription, the appropriate table will be created.
    If the class name has an underscore, the first part will be used as the group name.
    For example:
        Analysis_Dihedrals(IsDescription):
    will be a table 'dihedrals' in the 'analysis' group.
    """
    global _tables
    
    if os.path.exists(filename):
        print "Will not create datastore file at %s, file exists!" % filename
        return False
        
    for name in dir(module):
       obj = getattr(module, name)
       if inspect.isclass(obj) and obj.__base__.__name__ == 'IsDescription':
           _tables.append(obj)
    
    h5file = openFile(filename, mode="w", title=title)
    group = h5file.createGroup("/", 'detector', 'Detector information')
    table = h5file.createTable(group, 'readout', Particle, "Readout example")

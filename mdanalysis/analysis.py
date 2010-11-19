#!/usr/bin/python
# -*- coding: utf-8 -*-

from MDAnalysis import *
from MDAnalysis import collection, SelectionError

import tables
import numpy
import os

from nearby import NearbyListAnalysis

def main():
    from MDAnalysis.tests.datafiles import PSF,DCD
    print "Loading reference system: %s, %s" % (PSF, DCD)
    ref = Universe(PSF, DCD, permissive=True)
    print "Loading trajectory: %s" % (DCD)
    trj = Universe(PSF, DCD, permissive=True)
    
    # add the metadata for this analysis to the database
    analysis = Analysis('test.h5', readonly=False)
    
    # Test metadata
    analysis.add_metadata('/metadata/test', { 'a': '1', 'b': '2' })
    analysis.add_metadata('/metadata/test1', { 'a': '1', 'b': '2' })
    analysis.add_metadata('/metadata/test2', { 'a': '1', 'b': '2' })
    
    # Test timeseries
    analysis.add_timeseries('/timeseries/com/COM_ALL', Timeseries.CenterOfMass(ref.atoms))
    analysis.add_timeseries('/timeseries/com/COM_ALL1', Timeseries.CenterOfMass(ref.atoms))
    analysis.add_timeseries('/timeseries/com/COM_ALL2', Timeseries.CenterOfMass(ref.atoms))
    
    # Test sequence
    # analysis.add_to_sequence('/sequence/nearby_list', NearbyListAnalysis(), format=tables.Float32Atom(shape=()), array=True)
    
    analysis.run(trj=trj, ref=ref)
    analysis.save()
    analysis.close()

def split_path(path):
    split_path = path.split('/')
    leaf = split_path.pop()
    leaf_path = '/'.join(split_path)
    return (leaf_path, leaf)

class Analysis(object):
    _filename = None
    _readonly = True
    _title = None
    _h5f = None
    
    # nodes holds all Table and Array objects which take care of table creation and writing
    _nodes = {}
    
    # the following dicts store the actual analyses which are processed
    _sequential = {}
    _timeseries = {}
    
    def __init__(self, filename, title="datastore", readonly=True):
        self._filename = filename
        self._title = title
        self._readonly = readonly
        self._h5f = self.open_or_create()

    def get_or_create_array(self, path):
        if path not in self._nodes:
            self._nodes[path] = Array(self._h5f, path)
        return self._nodes[path]
            
    def get_or_create_column(self, path, format):
        (table_path, col_name) = split_path(path)
        if table_path not in self._nodes:
            self._nodes[table_path] = Table(self._h5f, table_path)
        return self._nodes[table_path].column(col_name, format)
    
    #analysis.add_metadata('/metadata/trajectory', { 'psf': psf_file, 'pdb': pdb_file, 'dcd': dcd_file, 'frames': num_frames, 'firsttimestep': first_timestep, 'dt': dt })
    def add_metadata(self, path, data, format=tables.StringCol(64)):
        #path is to the table
        #data has the columns
        
        print "Loading metadata..."
        for k, v in data.items():
            col_path = '%s/%s' % (path, k)
            col = self.get_or_create_column(col_path, format)
            col.load(str(v))
        print "Done."
    
    #analysis.add_timeseries('/protein/dihedrals/PEPA_139', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 139 N", "atom PEPA 139 CA", "atom PEPA 139 CB", "atom PEPA 139 CG")))
    def add_timeseries(self, path, timeseries, format=tables.Float32Col()):
        if path in self._timeseries:
            raise Exception('Timeseries with path %s already exists in this analysis!' % path)
        else:
            col = self.get_or_create_column(path, format)
            self._timeseries[path] = (timeseries, col)
   
    #analysis.add_to_sequence('/protein/rmsd/backbone', RMSD(ref, trj, selection='backbone')) 
    def add_to_sequence(self, path, processor, format=tables.Float32Col(), array=False):
        if path in self._sequential:
            raise Exception('Sequential processor with path %s already exists in this analysis!' % path)
        else:
            if array:
                node = self.get_or_create_array(path)
            else:
                node = self.get_or_create_column(path, format)
            self._sequential[path] = (processor, node)
    
    def run(self, trj, ref):
        self._trj = trj
        self._ref = ref
        
        if len(self._timeseries) > 0:
            print "Starting timeseries analysis..."
            collection.clear()
            for path, tpl in self._timeseries.items():
                print " Adding timeseries: %s" % path
                collection.addTimeseries(tpl[0])
            
            print " Computing..."
            collection.compute(self._trj.trajectory)
            print " Done computing."
        
            print "Loading data..."
            for i, path in enumerate(self._timeseries.keys()):
                print " loading table %s with %d values..." % (path, len(collection[i][0]))
                self._timeseries[path][1].load(list(collection[i][0]))
            print "Done timeseries analysis."
        
        if len(self._sequential) > 0:
            print "Running sequential analyses..."
            for path, tpl in self._sequential.items():
                print " Preparing %s" % path
                tpl[0].prepare(ref=self._ref, trj=self._trj)
            frames = self._trj.trajectory
            print " Processing %d frames..." % frames.numframes
            for i, f in enumerate(frames):
                if i % len(frames)/10 == 0:
                    print ".",
                for path, tpl in self._sequential.items():
                    tpl[0].process(f)
            print " done."
            print " Loading result data..."
            for path, tpl in self._sequential.items():
                tpl[1].load(tpl[0].results())
            print "Done sequential analysis."
        
    def save(self):
        print "Setting up and saving all tables and arrays..."
        for path, n in self._nodes.items():
            print " Node: %s" % path
            n.write()
        
    def close(self):
        print "Closing H5 file..."
        self._h5f.flush()
        self._h5f.close()
    
    def open_or_create(self):
        """ 
        Create the pytables data file at the given filename.
        """
        if not os.path.exists(self._filename):
            if self._readonly:
                raise Exception('Read-only open requested on file (%s) that doesn\'t exist!' % self._filename)
            mode = "w"
        else:
            # file exists
            if self._readonly:
                mode = "r"
            else:
                mode = "r+"
        return tables.openFile(self._filename, mode=mode, title=self._title)

class Column(object):
    _data = None # data to be written
    _dirty = False
    path = None
    name = None
    format = None
    
    def __init__(self, path, name, format):
        # print "Creating Column(%s)" % name
        self.path = path
        self.name = name
        self.format = format
        self._data = []
        self._dirty = False
    
    def load(self, data):
        # print "Loading data:"
        if type(data) is list:
            self._data += data
        else:
            self._data.append(data)
        self._dirty = True
    
    def dirty_row_count(self):
        return len(self._data)
    
    def next_dirty_row(self):
        if not self._dirty:
            return False
        
        row = self._data.pop(0)
        if len(self._data) == 0:
            self._dirty = False
        return row

class Array(Column):
    """ Array inherits from Column because it's basically a table with a single column """
    
    def __init__(self, h5f, full_path, format=tables.ObjectAtom()):
        self._h5f = h5f
        self.full_path = full_path
        (path, name) = split_path(self.full_path)
        super(Array, self).__init__(path, name, format)
    
    def setup(self):
        print "Setting up array at: %s" % (self.full_path)
        split_path = self.full_path.split('/')
        # traverse the path
        node = self._h5f.getNode('/')
        for n in split_path[1:]:
            parent = node._v_pathname
            path = (parent+'/'+n).replace('//','/')
            try:
                # try to get the node
                node = self._h5f.getNode(path)
            except tables.NoSuchNodeError:
                if path == self.full_path:
                    # we are at the array but it doesn't exist yet so create it
                    # node = self._h5f.createEArray(node, self.name, self.format, (len(self._data[0]), ), expectedrows=25000)
                    node = self._h5f.createVLArray(node, self.name, self.format, filters=tables.Filters(1))
                    print "Created array: %s" % path
                else:
                    # we are at a group that doesn't exist yet
                    node = self._h5f.createGroup(node, n)
                    print "Created group: %s" % path
        
        # node should now be the array node
        self._node = node
        return self._node
        
    def write(self):
        if self.dirty_row_count() == 0:
            print "Array %s has no rows to write, skipping it!" % self.full_path
            return False
        
        # first make sure the array is setup
        self.setup()
        
        print "Appending %d rows..." % self.dirty_row_count()
        while self._dirty:
            self._node.append(self.next_dirty_row())
        print " Done."

class Table(object):
    """ Dataset stored with PyTables
    """
    
    # PyTables table group, name and description
    
    _h5f = None # HDF5 file object (not initialized by this object)
    _node = None # the initialized table
    
    path = None
    name = None  # name of the table to be created/used
    
    _columns = {}
    
    def __init__(self, h5f, path):
        # print "Creating Table(%s)" % path
        self._h5f = h5f
        self._node = None
        self._columns = {}
        self.path = path
        self.name = self.path.split('/')[-1]
    
    def _description(self):
        desc = {}
        for col in self._columns.values():
            desc[col.name] = col.format
        return desc
    
    def column(self, name, format):
        # print "Adding column %s to table %s" % (name, self.path)
        if name not in self._columns:
            self._columns[name] = Column(self.path, name, format)
        return self._columns[name]
    
    def setup(self):
        print "Setting up table at: %s" % self.path
        split_path = self.path.split('/')
        # traverse the path
        node = self._h5f.getNode('/')
        for n in split_path[1:]:
            parent = node._v_pathname
            path = (parent+'/'+n).replace('//','/')
            try:
                # try to get the node
                node = self._h5f.getNode(path)
            except tables.NoSuchNodeError:
                if path == self.path:
                    # we are at the table but it doesn't exist yet so create it
                    node = self._h5f.createTable(node, n, self._description(), expectedrows=25000)
                    print "Created table: %s" % path
                else:
                    # we are at a group that doesn't exist yet
                    node = self._h5f.createGroup(node, n)
                    print "Created group: %s" % path
            else:
                # if we find the node but the descriptions differ (column(s) added or removed)
                if path == self.path and set(node.description._v_colObjects.keys()) != set(self._description().keys()):
                    print "Column(s) modified for table: %s" % self.path
                    print set(node.description._v_colObjects.keys()) ^ set(self._description().keys())
                    copy_node = self._h5f.createTable(node._v_parent, self.name+'_COPY', self._description(), expectedrows=25000)
                    node.attrs._f_copy(copy_node)
                    for i in xrange(node.nrows):
                        copy_node.row.append()
                    copy_node.flush()
                    # copy the data from the old table
                    for col in node.description._v_colObjects:
                        # only if the col is in the new table
                        if col in self._description():
                            getattr(copy_node.cols, col)[:] = getattr(node.cols, col)[:]
                    copy_node.flush()
                    node.remove()
                    copy_node.move(parent, self.name)
                    node = copy_node
                    print "Table %s updated with new column(s)." % self.path
        
        # node should now be the table object
        self._node = node
        return self._node
        
    def write(self):
        num_rows = [ col.dirty_row_count() for col in self._columns.values() ]
        num_rows = set(num_rows)
        if len(num_rows) > 1:
            raise Exception('Inconsistent number of rows to write: %s' % num_rows)
        num_rows = list(num_rows)[0]
        
        # first make sure the table is setup
        self.setup()
        
        print "Appending %d rows..." % num_rows
        row = self._node.row
        for i in range(num_rows):
            if i % 10 == 0:
                print '.',
            for col in self._columns.values():
                row[col.name] = col.next_dirty_row()
            row.append()
        print " Done."
        self._node.flush()
        
if __name__ == '__main__':
    main()

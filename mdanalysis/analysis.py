#!/usr/bin/python
# -*- coding: utf-8 -*-

from MDAnalysis import *
from MDAnalysis import collection, SelectionError

import tables
import numpy
import os

class Analysis(object):
    """
    
    """
    _filename = None
    _readonly = True
    _title = None
    _h5f = None
    
    # tables holds all Table objects which take care of table creation and writing
    _tables = {}
    
    # the following dicts store the actual analyses which are processed
    _sequential = {}
    _timeseries = {}
    
    def __init__(self, filename, title="datastore", readonly=True):
        self._filename = filename
        self._title = title
        self._readonly = readonly
        self._h5f = self.open_or_create()
            
    def get_column(self, path, format=tables.Float32Col()):
        split_path = path.split('/')
        col_name = split_path.pop()
        table_path = '/'.join(split_path)
        
        if table_path not in self._tables.keys():
            self._tables[table_path] = Table(self._h5f, table_path)
        
        # add the table if necessary and add the column
        return self._tables[table_path].column(col_name, format)
    
    #analysis.add_timeseries('/protein/dihedrals/PEPA_139', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 139 N", "atom PEPA 139 CA", "atom PEPA 139 CB", "atom PEPA 139 CG")))
    def add_timeseries(self, path, timeseries):
        if path in self._timeseries.keys():
            raise Exception('Timeseries with path %s already exists in this analysis!' % path)
        
        col = self.get_column(path)
        self._timeseries[path] = (timeseries, col)
   
    #analysis.add_to_sequence('/protein/rmsd/backbone', RMSD(ref, trj, selection='backbone')) 
    def add_to_sequence(self, path, processor):
        if path in self._sequential:
            raise Exception('Sequential processor with path %s already exists in this analysis!' % path)
        
        col = self.get_column(path)
        self._sequential[path] = (processor, col)
    
    def run(self, trj, ref):
        self._trj = trj
        self._ref = ref
        
        print "Starting timeseries analysis..."
        collection.clear()
        for path, tpl in self._timeseries.items():
            print " Adding timeseries: %s" % path
            collection.addTimeseries(tpl[0])
        print " Computing..."
        
        collection.compute(self._trj.dcd)
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
            ten_percent = int(float(frames.numframes)/10.0)
            for i, f in enumerate(frames):
                if i % ten_percent == 0:
                    print ".",
                for path, tpl in self._sequential.items():
                    tpl[0].process(f)
            print " done."
            print " Loading result data..."
            for path, tpl in self._sequential.items():
                tpl[1].load(tpl[0].results())
            print "Done sequential analysis."
        
    def save(self):
        print "Setting up and saving all tables..."
        for path, t in self._tables.items():
            print " Table: %s" % path
            t.setup()
            t.write()
    
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
            raise Exception('Tried to get row from a column without any data!')
        
        row = self._data.pop(0)
        if len(self._data) == 0:
            self._dirty = False
        return row
  
class Table(object):
    """ Dataset stored with PyTables
    """
    
    # PyTables table group, name and description
    
    _h5f = None # HDF5 file object (not initialized by this object)
    _table = None # the initialized table
    
    path = None
    name = None  # name of the table to be created/used
    
    _columns = {}
    
    def __init__(self, h5f, path):
        # print "Creating Table(%s)" % path
        self._h5f = h5f
        self.path = path
        self._table = None
        self._columns = {}
    
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
        split_path = self.path.split('/')
        table_name = split_path[-1]
        self.name = table_name
        
        path = []
        node = self._h5f.getNode('/')
        for n in split_path:
            try:
                # print "Looking for node: %s/%s" % ('/'.join(path), n)                
                node = self._h5f.getNode('/%s' % '/'.join(path), n)
                if n is table_name:
                    # loop through all columns.
                    #   if col exists, do nothing
                    #   if it doesn't, add it to the table
                    print "WARNING: adding new columns to tables not yet supported"
            except tables.NoSuchNodeError:
                if n is table_name:
                    # n is the table
                    # print "No table found, creating it..."
                    node = self._h5f.createTable(node, n, self._description(), expectedrows=25000)
                else:
                    # print "No group found, creating it..."
                    node = self._h5f.createGroup('/%s' % '/'.join(path), n)
            finally:
                path.append(n)
                if n is table_name:
                    self._table = node
        return self._table
        
    def write(self):
        num_rows = [ col.dirty_row_count() for col in self._columns.values() ]
        num_rows = set(num_rows)
        if len(num_rows) > 1:
            raise Exception('Inconsistent number of rows to write: %s' % num_rows)
        num_rows = list(num_rows)[0]
        
        print "Appending %d rows..." % num_rows
        row = self._table.row
        ten_percent = int(float(num_rows)/10.0)
        for i in range(num_rows):
            if i % ten_percent == 0:
                print ".",
            for col in self._columns.values():
                row[col.name] = col.next_dirty_row()
            row.append()
        print " Done."
        self._table.flush()
    
    # def __repr__(self):
    #        return '<'+self.__class__.__name__+' '+repr(self.path)+'/'+repr(self.name)+'>'
    
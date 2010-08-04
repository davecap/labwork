#!/usr/bin/python
# -*- coding: utf-8 -*-

from tables import *
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
    #   'path' => (class, post_processor)
    _sequential = {}
    #   'path' => (class, post_processor)
    _timeseries = {}
    
    def __init__(self, filename, title="datastore", readonly=True):
        self._filename = filename
        self._title = title
        self._readonly = readonly
        self._h5f = self.open_or_create()
            
    # get the table, "create" it if necessary
    #   create means just create the Table class object, won't create the table until run()
    def table(path):
        if path not in self._tables:
            self._tables[path] = Table(self._h5f, path)
        return self._tables[path]
    
    def column(path, format=Int32col()):
        split_path = path.split('/')
        col = split_path.pop()
        group_table = '/'.join(split_path)

        # add the table if necessary and add the column
        return self.table(group_table).column(col, format)
    
    #analysis.add_timeseries('/protein/dihedrals/PEPA_139', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 139 N", "atom PEPA 139 CA", "atom PEPA 139 CB", "atom PEPA 139 CG")), pp=(lambda x: x*180./pi))
    def add_timeseries(path, timeseries, pp=None):
        if path in self._timeseries:
            raise Exception('Timeseries with path %s already exists in this analysis!' % path)
        
        col = self.column(path)
        self._timeseries[path] = (timeseries, col)
   
    #analysis.add_to_sequence('/protein/rmsd/backbone', RMSD(ref, trj, selection='backbone')) 
    def add_to_sequence(path, processor):
        if path in self._sequential:
            raise Exception('Sequential processor with path %s already exists in this analysis!' % path)
        
        col = self.column(path)
        self._sequential[path] = (processor, col)
    
    def run(self, trj):
        print "Running sequential analyses..."
        for path, tpl in self._sequential.items():
            print "Preparing %s" % path
            tpl[0].prepare()
        
        frames = trj.trajectory
        for ts in frames:
            for a in analyses:
                a.process(ts)
        for path, tpl in self._sequential.items():
            tpl[1].load(tpl[0].results())
        
        # open the datastore
         for ts in self._timeseries:
             # validate datastore against all analyses to be done
             # create dataset if necessary
             collection.addTimeseries(ts[0])
        
         #data = universe.dcd.correl(collection, stop=5)
         collection.compute(trj.dcd)
        
         for i, ts in enumerate(self._timeseries):
             # write to the datastore
             ts[1].load(collection[i][0])
        
    def save(self):
        #save all the tables/columns
        pass
    
    def close(self):
        pass
    
    def write(self):
        for t in self._tables:
            t.write()
    
    def close(self):
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
    _data = None # data to ve written
    _rows = 0 # rows to be written
    _dirty = False
    path = None
    name = None
    format = None
    
    def __init__(self, path, name, format):
        self.path = path
        self.name = name
        self.format = format
        self._data = []
        self._dirty = False
    
    def load(self, data):
        self._data.append(data)
        self._dirty = True
    
    def save(self):
        self._dirty = False
  
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
        if not h5f or not h5f.isopen():
            raise Exception('Closed H5 file descriptor passed to Analysis class: %s' % self)
        self._h5f = h5f
        self.path = path
        else:
            raise Exception('Too many levels in the table path: %s' % path)
    
    def _description(self):
        desc = {}
        for col in self._columns.values():
            desc[col.name] = col.format
        return desc
    
    def column(self, name, format):
        # get column, make the object if necessary
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
                print "Looking for node: /%s/%s" % ('/'.join(path), n)                
                node = self._h5f.getNode('/%s' % '/'.join(path), n)
                if n is table_name:
                    # loop through all columns.
                    #   if col exists, do nothing
                    #   if it doesn't, add it to the table
                    print "WARNING: adding new columns to tables not yet supported"
            except tables.NoSuchNodeError:
                if n is table_name:
                    # n is the table
                    print "No table found, creating it..."
                    node = self._h5f.createTable(node, t, self._description(), expectedrows=100)
                else:
                    print "No group found, creating it..."
                    node = self._h5f.createGroup('/%s' % '/'.join(path), n)
            finally:
                path.append(n)
                if n is table_name:
                    self._table = node
        return self._table
        
    def write(self, data, col=None):
        """
        Write the data in _data to the table.
        _data should be a dict with the format:
            _data['<table column name>'] = numpy.array([...])
        """
        row = self._table.row
        if col:
            for i in range(len(self._data[col])):
                row[col] = self._data[col][i]
                row.append()
        else:
            for i in range(self._rows):
                for col in self._data.keys()
                    row[col] = self._data[col][i]
                row.append()
        self._table.flush()
    
    def __repr__(self):
        return '<'+self.__class__.__name__+' '+repr(self.name)+', '+repr(self.id)+', '+repr(self._table_group)+'/'+repr(self._table_name)+'>'

class SequentialAnalysis(object):
    """
    Performs analysis per frame in a trajectory.
    Puts the results in the _data dictionary.
    """
    
    def process(self, frame):
        raise NotImplementedError()

    def results(self):
        raise NotImplementedError()
    
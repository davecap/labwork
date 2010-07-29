#!/usr/bin/python
# -*- coding: utf-8 -*-

from tables import *
import numpy
import os

class Analysis(object):
    """
        
    """

class Analysis(object):
    """
    
    """
    _filename = None
    _readonly = True
    _title = None
    _h5f = None
    
    _tables = {}
    
    def __init__(self, filename, title="datastore", readonly=True):
        self._filename = filename
        self._title = title
        self._readonly = readonly
        
        self._h5f = self.open_or_create()
    
    def _add_column(path):
        # get or create the table
        #   add the column to the description
    
    def add_timeseries(path, timeseries, pp=None):
        #analysis.add_timeseries('/protein/dihedrals/PEPA_139', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 139 N", "atom PEPA 139 CA", "atom PEPA 139 CB", "atom PEPA 139 CG")), pp=(lambda x: x*180./pi))
        split_path = path.split('/')
 
    
    # def table(self, group, name):
    #     """
    #     Add a table to the DataStore. Create it if necessary and prepare for writing.
    #     Pass this function a Tabular object.
    #     """
    #     path = '/%s/%s' % (group, name)
    #     try:
    #         return self._tables[path]
    #     except Exception:
    #         # table not found
    #         pass
    #     
    #     self._tables[path] = new_table
    #     return new_table
    
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
        
class DataStoreTable(object):
    """ Dataset stored with PyTables
    """
    
    # PyTables table group, name and description
    
    _h5f = None # HDF5 file object (not initialized by this object)
    _table = None # the initialized table
    
    name = None  # name of the table to be created/used
    group = None # group of this table, will be created if it doesn't exist
    description = { }
    # Example:
    #   description = {
    #       'time': Int32Col(),
    #       'name': StringCol(64),
    #       'count': Int32Col(),
    #   }
    
    _data = { } # data to be written
    _rows = 0
    
    def __init__(self, h5f):
        if not h5f or not h5f.isopen():
            raise Exception('Closed H5 file descriptor passed to Analysis class: %s' % self)
        self._h5f = h5f
        
        try:
            # look for the group
            group_node = self._h5f.getNode('/', self.group)
        except tables.NoSuchNodeError:
            print "No group found for /%s, creating it..." % (self.group)
            group_node = self._h5f.createGroup('/', self.group)
        
        try:
            # look for the table
            self._table = self._h5f.getNode('/%s' % self.group, self.name)
        except tables.NoSuchNodeError:
            print "No table found for /%s/%s, creating it..." % (self.group, self.name)
            self._table = self._h5f.createTable(group_node, self.name, self.description, self.name)
        
    def write(self, col=None):
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
    
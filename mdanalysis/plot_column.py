#!/usr/bin/python
# -*- coding: utf-8 -*-

import optparse
import numpy
import numpy.linalg
import scipy.stats
from MDAnalysis import *
from MDAnalysis import collection, SelectionError
import MDAnalysis.core.rms_fitting
from MDAnalysis.core.AtomGroup import Residue, AtomGroup
from math import pi

from analysis import *

import tables

def main():
    usage = """
        usage: %prog [options] <H5 File> <path>
    """
    parser = optparse.OptionParser(usage)
    
    options, args = parser.parse_args()
    
    if len(args) < 2:
        parser.error("No input file specified or path specified")
    
    h5_file = args[0]
    path = args[1]
    
    column = path.split('/')[-1]
    ps = path.split('/')
    ps.pop()
    table = '/'.join(ps)
    
    import matplotlib.pyplot as plt
    
    h5f = tables.openFile(h5_file, mode="r")
    tbl = h5f.getNode(table)
    data = tbl.read(field=column)
    h5f.close()
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(data)
    ax.set_xlabel(r'Frame')
    ax.set_ylabel(column)
    ax.set_title(path)
    plt.show()
    

if __name__ == '__main__':
    main()

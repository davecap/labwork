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
import matplotlib.pyplot as plt

from analysis import *

import tables

def main():
    usage = """
        usage: %prog [options] <H5 File> <path1> <path2> ... <pathN>
    """
    parser = optparse.OptionParser(usage)
    
    options, args = parser.parse_args()
    
    try:
        h5_file = args[0]
    except:
        parser.error("No input file specified")
    
    del args[0]
    if len(args) == 0:
        print "No path specified, showing all possible paths..."
        h5f = tables.openFile(h5_file, mode="r")
        for g in h5f.root._v_children.keys():
            for t in h5f.root._v_children[g]._v_children.keys():
                for c in h5f.root._v_children[g]._v_children[t].description._v_names:
                    print "/%s/%s/%s" % (g, t, c)
        h5f.close()
        parser.error("No path specified")
        
    h5f = tables.openFile(h5_file, mode="r")
    fig = plt.figure()
    ax = fig.add_subplot(111)
        
    for path in args:
        column = path.split('/')[-1]
        ps = path.split('/')
        ps.pop()
        table = '/'.join(ps)
        tbl = h5f.getNode(table)
        ax.plot(tbl.read(field=column), label=path)
        
    h5f.close()
    
    ax.legend()
    ax.set_xlabel(r'Frame')
    ax.set_ylabel(r'')
    ax.set_title(h5_file)
    plt.show()
    

if __name__ == '__main__':
    main()

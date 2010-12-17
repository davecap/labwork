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
    
    try:
        h5_file = args[0]
    except:
        parser.error("No input file specified")
    
    try:
        path = args[1]
    except:
        print "No path specified, showing all possible paths..."
        h5f = tables.openFile(h5_file, mode="r")
        for g in h5f.root._v_children.keys():
            for t in h5f.root._v_children[g]._v_children.keys():
                try:
                    for c in h5f.root._v_children[g]._v_children[t].description._v_names:
                        print "/%s/%s/%s" % (g, t, c)
                except:
                    print "/%s/%s" % (g, t)
        h5f.close()
        parser.error("No path specified")
    
    column = path.split('/')[-1]
    ps = path.split('/')
    ps.pop()
    table = '/'.join(ps)
    
    import matplotlib.pyplot as plt
    
    h5f = tables.openFile(h5_file, mode="r")
    tbl = h5f.getNode(table)
    data = tbl.read(field=column)
    h5f.close()

    rad2deg = (lambda x: x*180./pi)
    data_deg = []
    for rad in data:
        deg = rad2deg(rad)
        if deg < 0:
            deg += 360
        data_deg.append(deg)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(data_deg)
    ax.set_xlabel(r'Frame')
    ax.set_ylabel(column)
    ax.set_title(path)
    plt.show()
    

if __name__ == '__main__':
    main()

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
                for c in h5f.root._v_children[g]._v_children[t].description._v_names:
                    print "/%s/%s/%s" % (g, t, c)
        h5f.close()
        parser.error("No path specified")
    
    import matplotlib.pyplot as plt
    
    h5f = tables.openFile(h5_file, mode="r")
    arr = h5f.getNode(path)
    rows = arr.read()
    h5f.close()
    
    # for each row: [ ('name', position ) ... ]
    positions = {}
    for r in rows:
        # put it in a dict: positions['name'] = []
        for d in r:
            n = d[0].split(':')[0]
            if n in positions:
                positions[n].append(float(d[1]))
            else:
                positions[n] = [float(d[1])]
    
    # plot histogram for each position
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for n, p in positions.items():
        ax.hist(p, bins=200, label=n)
    ax.legend()
    ax.set_xlabel(r'Distance along cylinder')
    ax.set_ylabel('Frequency')
    ax.set_title(path)
    plt.show()
    

if __name__ == '__main__':
    main()

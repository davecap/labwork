#!/usr/bin/python
# -*- coding: utf-8 -*-

import optparse
import numpy
from scipy.stats import sem

from math import pi
from configobj import ConfigObj, flatten_errors
import sys
import tables

import os
import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as plt
import matplotlib.font_manager

def rad2deg(rad):
    deg = rad*180./pi
    if deg < 0:
        deg += 360.
    return deg
    
def get_all_fields(h5_file):
    fields = []

    return fields

def main():
    usage = """
        usage: %prog [options] <config.ini files>
    """
    parser = optparse.OptionParser(usage)
    parser.add_option("-x", dest="x_column", default=None, help="X Column REQUIRED")    
    # parser.add_option("-y", dest="y_column", default=None, help="Y Column REQUIRED")
    options, args = parser.parse_args()    
    
    _h5_filename = 'analysis.h5'
    
    if not len(args):
        parser.error("No config files specified!")
    elif not options.x_column:
        # get the first config
        # show all paths
        config = ConfigObj(args[0])
        project_path = os.path.abspath(os.path.dirname(config.filename))
        print "A path is required! showing all possible paths..."
        for r_id, r_config in config['replicas'].items():
            try:
                h5_file = os.path.join(project_path, r_id, _h5_filename)
                h5f = tables.openFile(h5_file, mode="r")
                for g in h5f.root._v_children.keys():
                    for t in h5f.root._v_children[g]._v_children.keys():
                        for c in h5f.root._v_children[g]._v_children[t].description._v_names:
                            print "/%s/%s/%s" % (g, t, c)
                h5f.close()
                break
            except:
                continue
        parser.error("X path is required.")
    
    # ok to continue...
    h5_field = options.x_column.split('/')[-1]
    ps1 = options.x_column.split('/')
    ps1.pop()
    h5_table = '/'.join(ps1)
    vfunc = numpy.vectorize(rad2deg)
    
    # set up the plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'Chi1 of Residue 139')
    ax.set_ylabel(ps1[-1])
    ax.set_title(options.x_column)
    
    for config_file in args:
        if not os.path.exists(config_file):
            sys.stderr.write("Config file not found at: %s" % config_file)
            continue
    
        config = ConfigObj(config_file)
        project_path = os.path.abspath(os.path.dirname(config.filename))
        data = { 'x': numpy.array([]), 'y': numpy.array([]), 'yerr': numpy.array([]) }
    
        # loop through replicas in the config file
        for r_id, r_config in config['replicas'].items():
            try:
                h5_file = os.path.join(project_path, r_id, _h5_filename)
                h5f = tables.openFile(h5_file, mode="r")
                tbl = h5f.getNode(h5_table)
                field_data = tbl.read(field=h5_field)
                if 'dihedral' in options.x_column:
                    field_data = vfunc(field_data)
            except Exception, e:
                sys.stderr.write("Error reading H5 file (%s): %s\n" % (h5_file,e))
            else:
                data['y'] = numpy.append(data['y'], numpy.mean(field_data))
                data['x'] = numpy.append(data['x'], float(r_config['coordinate']))
                data['yerr'] = numpy.append(data['yerr'], sem(field_data))
            finally:
                h5f.close()
    
        # plot the current dataset
        ax.errorbar(data['x'], data['y'], yerr=data['yerr'], fmt='o', label=config['title'])
    ax.legend()
    plt.savefig(h5_field + '.eps')
    

if __name__ == '__main__':
    main()

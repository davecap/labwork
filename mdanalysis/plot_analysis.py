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

def main():
    usage = """
        usage: %prog [options]
    """
    parser = optparse.OptionParser(usage)
    parser.add_option("-x", dest="x_column", default=None, help="X Column REQUIRED")    
    # parser.add_option("-y", dest="y_column", default=None, help="Y Column REQUIRED")
    options, args = parser.parse_args()    
    
    _h5_filename = 'analysis.h5'
        
    if options.x_column:
        h5_field = options.x_column.split('/')[-1]
        ps1 = options.x_column.split('/')
        ps1.pop()
        h5_table = '/'.join(ps1)
        vfunc = numpy.vectorize(rad2deg)
        
        # set up the plot
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(r'Replica')
        ax.set_ylabel(options.x_column)
    
    for config_file in args:
        if not os.path.exists(config_file):
            sys.stderr.write("Config file not found at: %s" % config_file)
            continue
    
        config = ConfigObj(config_file)
        project_path = os.path.abspath(os.path.dirname(config.filename))
        data = { 'x': numpy.array([]), 'y': numpy.array([]), 'yerr': numpy.array([]) }
    
        # loop through replicas in the config file
        for r_id, r_config in config['replicas'].items():
            # check if h5 file exists for the current replica
            h5_file = os.path.join(project_path, r_id, _h5_filename)
            if not os.path.exists(h5_file):
                #sys.stderr.write("H5 file not found... skipping %s/%s\n" % (config['title'], r_id))
                continue
            else:
                #sys.stderr.write("Reading data from: %s/%s\n" % (config['title'], r_id))
                if not options.x_column:        
                    print "A path is required! showing all possible paths..."
                    h5f = tables.openFile(h5_file, mode="r")
                    for g in h5f.root._v_children.keys():
                        for t in h5f.root._v_children[g]._v_children.keys():
                            for c in h5f.root._v_children[g]._v_children[t].description._v_names:
                                print "/%s/%s/%s" % (g, t, c)

                    h5f.close()
                    parser.error("X path is required.")
                else:
                    # read the actual data
                    h5f = tables.openFile(h5_file, mode="r")
                    try:
                        tbl = h5f.getNode(h5_table)
                        field_data = tbl.read(field=h5_field)
                        if 'dihedral' in options.x_column:
                            field_data = vfunc(field_data)
                    except Exception, e:
                        sys.stderr.write("Error reading analysis file... %s\n" % e)
                        # data = numpy.append(data, 0.)
                    else:
                        data['y'] = numpy.append(data['y'], numpy.mean(field_data))
                        data['x'] = numpy.append(data['x'], float(r_config['coordinate']))
                        data['yerr'] = numpy.append(data['yerr'], sem(field_data))
                    finally:
                        h5f.close()
    
        # plot the current dataset
        ax.errorbar(data['x'], data['y'], yerr=data['yerr'], fmt='o', label=config['title'])
    # ax.set_title(path)
    # plt.show()
    ax.legend()
    plt.savefig(h5_field + '.eps')

    

if __name__ == '__main__':
    main()

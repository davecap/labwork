#!/usr/bin/python
# -*- coding: utf-8 -*-

import optparse
import numpy
import os
import sys
import tables
from configobj import ConfigObj, flatten_errors

import matplotlib
import matplotlib.pyplot as plt

def extract(config, path):
    project_path = os.path.abspath(os.path.dirname(config.filename))
    
    # loop through replicas in the config file
    for r_id, r_config in config['replicas'].items():
        #print "Extracing data from replica %s" % str(r_id)
        try:
            h5_file = os.path.join(project_path, r_id, 'analysis.h5')
            h5f = tables.openFile(h5_file, mode="r")
            arr = h5f.getNode(path)
            rows = arr.read()
            h5f.close()
        except Exception, e:
            sys.stderr.write("Error reading H5 file (%s): %s\n" % (h5_file,e))
            yield (r_config['coordinate'], [])
        else:
            yield (r_config['coordinate'], rows)
            
# try:
#     path = args[1]
# except:
#     print "No path specified, showing all possible paths..."
#     h5f = tables.openFile(h5_file, mode="r")
#     for g in h5f.root._v_children.keys():
#         for t in h5f.root._v_children[g]._v_children.keys():
#             for c in h5f.root._v_children[g]._v_children[t].description._v_names:
#                 print "/%s/%s/%s" % (g, t, c)
#     h5f.close()
#     parser.error("No path specified")
            

def get_positions(rows):
    positions = {}
    # for each row: [ ('name', position ) ... ]
    for r in rows:
        # put it in a dict: positions['name'] = []
        for d in r:
            n = d[0].split(':')[0]
            if n in positions:
                positions[n].append(float(d[1]))
            else:
                positions[n] = [float(d[1])]
    return positions

def water_wire_length(row, cutoff=5.0):
    sorted_row = sorted([ float(r[1]) for r in row if 'TIP3' in r[0] and float(r[1]) > 0 ])
    wire = [sorted_row[0]]
    for r in sorted_row[1:]:
        if r-wire[-1] <= cutoff:
            wire.append(r)
        else:
            break
    return wire[-1] - wire[0]
    
def main():
    usage = """
        usage: %prog [options] <config.ini files>
    """
    parser = optparse.OptionParser(usage)
    parser.add_option("-x", dest="x_column", default=None, help="X Column REQUIRED")    
    options, args = parser.parse_args()    
        
    if not len(args):
        parser.error("No config files specified!")
    
    #
    # CoM Differences for TIP3
    #
    
    #differences = []
    # for arg in args:
    #     print "Extracting data from config %s" % arg
    #     for rows in extract(arg, options.x_column):
    #         for r in rows:
    #             positions = [ float(i[1]) for i in r if 'TIP3' in i[0] ]
    #             positions.sort()
    #             prev = positions[0]
    #             for p in positions[1:]:
    #                 differences.append(p-prev)
    #                 prev = p
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.set_xlim(0,12)
    # 
    # (hist, bin_edges) = numpy.histogram(differences, bins=100, normed=True)
    # ax.plot(bin_edges[:-1], hist)
    # ax.set_xlabel(r'Distance between center of mass')
    # ax.set_ylabel('Frequency')
    # ax.set_title('TIP3 CoM Distance Distribution')
    #plt.show()
    
    # 
    # Water wire length distribution
    #
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    cmap = matplotlib.cm.get_cmap('Blues')
    pot_cmap = matplotlib.cm.get_cmap('OrRd')
    
    # for each system
    for arg in args:
        print "Extracting data from config %s" % arg
        if not os.path.exists(arg):
            sys.stderr.write("Config file not found at: %s" % arg)
            continue
        config = ConfigObj(arg)
        lengths = []
        # extract data from data files for each replica (coordinate)
        for (coord,rows) in extract(config, options.x_column):
            positions = []
            pot_positions = []
                
            # TIP3 + POT distributions
            
            for r in rows:
                positions += [ float(i[1]) for i in r if 'TIP3' in i[0] and float(i[1]) >= 0 ]         
                pot_positions += [ float(i[1]) for i in r if 'POT' in i[0] and float(i[1]) >= 0 ]                
            try:
                (hist, bin_edges) = numpy.histogram(positions, bins=100, normed=True)
                ax.scatter([ float(coord) for i in hist ], bin_edges[:-1], c=hist, cmap=cmap, s=20, marker='o', edgecolors="none")
                #ax.plot(bin_edges[:-1], hist, label=coord)
            except:
                 print "No data to plot for %s" % arg
            try:
                (hist, bin_edges) = numpy.histogram(pot_positions, bins=100, normed=True)
                ax.scatter([ float(coord) for i in hist ], bin_edges[:-1], c=hist, s=20, marker='o', cmap=pot_cmap, alpha=0.15, edgecolors="none")
                # ax.plot(bin_edges[:-1], hist, label=config['title'])
            except:
                print "No data to plot for %s" % arg
            
            
            #if len(rows) > 0:
            #    lengths += map(water_wire_length, rows)
            #    # (hist, bin_edges) = numpy.histogram(map(water_wire_length, rows), bins=30, normed=True)
        #(hist, bin_edges) = numpy.histogram(lengths, bins=50, normed=True)
        #ax.plot(bin_edges[:-1], hist, label=config['title'])

    # ax.legend()
    
    ax.set_xlabel(r'Chi1 of 139')
    # ax.set_ylabel('Water wire length distribution')
    ax.set_title(config['title']+' TIP3/POT distribution along cylinder vs Chi1 of residue 139')
    
    #ax.set_xlabel(r'Distance along cylinder for TIP3 ' + options.x_column)
    # ax.set_ylabel('Probability density of potassium')
    # ax.set_ylabel(r'Probability of > 15 A water wire from 132 to 286')    
    ax.set_ylim(0, 35)
    ax.set_xlim(0,360)
    plt.show()
    
    #
    # plot histogram for each position
    #
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111)

    # for n, p in positions.items():
    #     (hist, bin_edges) = numpy.histogram(p, bins=80, normed=True)
    #     # bins = ax.hist(p, bins=60, label=n)
    #     ax.plot(bin_edges[:-1], hist, label=n)
    # ax.legend()
    # ax.set_xlabel(r'Distance along cylinder')
    # ax.set_ylabel('Frequency')
    # ax.set_title(h5_file+path)
    # fname = (h5_file+path).replace('./','').replace('/','_')
    # plt.savefig('arrays/'+fname+'.eps')

if __name__ == '__main__':
    main()

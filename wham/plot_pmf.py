#!/usr/bin/python

import os
import optparse
from configobj import ConfigObj, flatten_errors
import tables
from math import pi
import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as plt
import matplotlib.font_manager
import numpy

def main():    
    usage = """
        usage: %prog [options] <key>,<wham outfile> <key2>,<wham outfile2>
    """
    
    parser = optparse.OptionParser(usage)
    # parser.add_option("-c", "--config", dest="config_file", default="config.ini", help="Config file [default: %default]")
    (options, args) = parser.parse_args()

    # if not os.path.exists(options.config_file):
    #     raise Exception('No config file found!')
    # config = ConfigObj(options.config_file)
    # project_path = os.path.abspath(os.path.dirname(config.filename))
    
    # setup the plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for f in args:
        print "Parsing arg %s..." % f
        if ',' in f:
            split_f = f.split(',')
            title = split_f[0]
            filename = split_f[1]
            #color = split_f[2]
            print "\tfound title: %s" % title
            #print "\tfound color: %s" % color
        else:
            title = f
            filename = f
        
        x = []
        y = []
        e = []

        infile = open(filename, 'r')
        while infile:
            line = infile.readline().strip()
            if len(line) == 0:
                break
            elif line.startswith('#'):
                continue
            split_line = [ s.strip() for s in line.split('\t') ]
            if len(split_line) > 3 and split_line[1] != 'inf':
                x.append(float(split_line[0]))
                y.append(float(split_line[1]))
                if split_line[2] == "nan":
                    e.append(0.0)
                else:
                    e.append(float(split_line[2]))
            else:
                print "Skipping line: %s" % line
        plt.errorbar(x, y, yerr=e, label=title)
        
    # complete and show the plot
    ax.set_xlim(0, 360)
    ax.set_ylim(0, 20)
    prop = matplotlib.font_manager.FontProperties(size=8)
    ax.legend(prop=prop)
    ax.set_xlabel(r'Coordinate')
    ax.set_ylabel(r'dG')
    ax.set_title('PMF of residue 139 dihedral angle for three CcO variants')
    ax.grid(True)
    plt.savefig('pmf.ps')
    #plt.show()
    
if __name__=='__main__':
    main()

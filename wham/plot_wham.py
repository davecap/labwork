#!/usr/bin/python

import os
import optparse
from configobj import ConfigObj, flatten_errors
import tables
from math import pi
import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as plt

import numpy

#from pydr import setup_config

def main():    
    usage = """
        usage: %prog [options]
    """
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-c", "--config", dest="config_file", default="config.ini", help="Config file [default: %default]")
    (options, args) = parser.parse_args()

    if not os.path.exists(options.config_file):
        raise Exception('No config file found!')
    
    #config = setup_config(options.config_file, create=False)
    config = ConfigObj(options.config_file)
    project_path = os.path.abspath(os.path.dirname(config.filename))

    wham_metadata_filename = os.path.join(project_path, 'wham_metadata')
    if os.path.exists(wham_metadata_filename):
        os.rename(wham_metadata_filename, wham_metadata_filename+'.backup')
    wham_metadata_file = open(wham_metadata_filename, 'w')
    
    # setup the plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for r_id, r_config in config['replicas'].items():
        print "Extracting data from replica: %s" % r_id
        h5_file = os.path.join(project_path, r_id, 'analysis.h5')
        if not os.path.exists(h5_file):
            print "H5 file not found... skipping this replica"
            continue
        coordinate = float(r_config['coordinate'])
        k = float(r_config['k'])
        
        h5f = tables.openFile(h5_file, mode="r")
        try:
            tbl = h5f.getNode('/protein/dihedrals')
        except Exception, e:
            print "Error reading analysis file... %s" % e
            h5f.close()
            continue
        
        data_chi1 = tbl.read(field='PEPA_139_CHI1')
        #data_chi2 = tbl.read(field='PEPA_139_CHI2')
        h5f.close()

        #rad2deg = (lambda x: x*1)
        rad2deg = (lambda x: x*180./pi)

        data_file_path = os.path.join(project_path, r_id, 'wham_data')
        if os.path.exists(data_file_path):
            os.rename(data_file_path, data_file_path+'.backup')
        f = open(data_file_path, 'w')
        degrees = []
        for rad in data_chi1:
            deg = rad2deg(rad)
            if deg < 0:
                deg += 360
            # write to file: "0.0    deg"
            f.write('0.0    %f\n' % deg)
            degrees.append(deg)
        f.close()
        wham_metadata_file.write('%s %s %s\n' % (data_file_path, coordinate, k))
        
        # generate umbrellas for each coord/k
        x = numpy.arange(coordinate-10.0,coordinate+10.1,0.1)
        y = [ 0.5*k*(i-coordinate)*(i-coordinate) for i in x ]
        plt.plot(x,y)
        # add degree data for each coord
        #y = [ 0.5*k*(i-coordinate)*(i-coordinate) for i in degrees ]
        #plt.scatter(degrees, y)
        
    wham_metadata_file.close()
    
    # complete and show the plot
    ax.set_xlim(0, 360)
    ax.set_ylim(0, 5)
    ax.set_xlabel(r'Coordinate')
    ax.set_ylabel(r'Harmonic Restraint')
    ax.set_title('')
    # ax.grid(True)
    plot_path = os.path.join(project_path, 'wham_plot.ps')
    plt.savefig(plot_path)
    #plt.show()
    
if __name__=='__main__':
    main()

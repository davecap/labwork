#!/usr/bin/python

import os
import sys
import optparse
from configobj import ConfigObj, flatten_errors
import tables
from math import pi
import tmpfile
import random

def main():    
    usage = """
        usage: %prog [options]
    """
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-c", "--config", dest="config_file", default="config.ini", help="Config file [default: %default]")
    parser.add_option("-p", "--percent", dest="percent", default=100, type="int", help="Percent of data to extract for WHAM [default: %default %%]")
    parser.add_option("-r", "--random", dest="random", default=False, action="store_true", help="Select data points at random [default: %default]")
    parser.add_option("-o", "--output-dir", dest="output_dir", default=None, help="Output directory [default: temporary dir]")
    parser.add_option("-s", "--start-index", dest="start_index", default=0, type="int", help="Starting index for column data [default: %default]")
    
    (options, args) = parser.parse_args()

    _wham_metadata_filename = 'wham_metadata'
    _h5_filename = 'analysis.h5'
    _h5_table = '/protein/dihedral'
    _h5_field = 'PEPA_139_CHI1'

    if not os.path.exists(options.config_file):
        raise Exception('No config file found!')
    
    config = ConfigObj(options.config_file)
    project_path = os.path.abspath(os.path.dirname(config.filename))
    
    if options.output_dir is None:
        # use a temp dir
        output_dir = tmpfile.mkdtemp()
        sys.stderr.write("Creating temporary dir at %s" % output_dir)
    elif options.output_dir.startswith('/'):
        # absolute path
        output_dir = options.output_dir
    else:
        # relative path
        output_dir = os.path.join(project_path, options.output_dir)
    
    # open the wham metadata file for writing
    wham_metadata_file = open(os.path.join(output_dir, _wham_metadata_filename), 'w')
    
    wham_metadata_file.write('# WHAM Metadata File generated from Python')
    wham_metadata_file.write('# Project Path: %s' % project_path)
    wham_metadata_file.write('# Config: %s' % options.config)
    wham_metadata_file.write('# Output Dir: %s' % output_dir)
    wham_metadata_file.write('# Percent Data Used: %s%' % str(options.percent))
    wham_metadata_file.write('# Random Selection: %s' % str(options.random))
    
    # loop through replicas in the config file
    for r_id, r_config in config['replicas'].items():
        sys.stderr.write("Extracting data from replica: %s" % r_id)
        
        # check if h5 file exists for the current replica
        h5_file = os.path.join(project_path, r_id, _h5_filename)
        if not os.path.exists(h5_file):
            sys.stderr.write("H5 file not found... skipping this replica")
            continue
        
        h5f = tables.openFile(h5_file, mode="r")
        try:
            tbl = h5f.getNode(_h5_table)
        except Exception, e:
            sys.stderr.write("Error reading analysis file... %s" % e)
            h5f.close()
            continue
        
        field_data = tbl.read(field=_h5_field)
        h5f.close()
        
        # extract a subset of the data
        sample_size = int(round(float(options.percent)/100.0 * float(len(field_data))))
        
        if options.random:
            field_data = random.sample(field_data[options.start_index:], sample_size)
        else:
            # slice data from start_index to sample_size
            if sample_size > len(field_data)-options.start_index:
                raise Exception('Required sample size (%d) is larger than dataset (%d) starting from index %d' % (sample_size, len(field_data), options.start_index))
            field_data = field_data[options.start_index:options.start_index+sample_size]
        
        #rad2deg = (lambda x: x*1)
        rad2deg = (lambda x: x*180./pi)

        # data_file_path = os.path.join(project_path, r_id, 'wham_data')
        data_file_path = os.path.join(output_dir, 'wham_data_%s' % r_id)
        if os.path.exists(data_file_path):
            os.rename(data_file_path, data_file_path+'.backup')
        
        f = open(data_file_path, 'w')
        for rad in field_data:
            deg = rad2deg(rad)
            if deg < 0:
                deg += 360
            f.write('0.0    %f\n' % deg)
        f.close()
        
        wham_metadata_file.write('%s %s %s\n' % (data_file_path, r_config['coordinate'], r_config['k']))
    wham_metadata_file.close()
    
if __name__=='__main__':
    main()

#!/usr/bin/python

import os
import optparse
from configobj import ConfigObj, flatten_errors
import tables
from math import pi

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
    for r_id, r_config in config['replicas'].items():
        print "Extracting data from replica: %s" % r_id
        h5_file = os.path.join(project_path, r_id, 'analysis.h5')
        if not os.path.exists(h5_file):
            print "H5 file not found... skipping this replica"
            continue
        r_config['coordinate']
        r_config['k']
        
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
        for rad in data_chi1:
            deg = rad2deg(rad)
            if deg < 0:
                deg += 360
            # write to file: "0.0    deg"
            f.write('0.0    %f\n' % deg)
        f.close()
        wham_metadata_file.write('%s %s %s\n' % (data_file_path, r_config['coordinate'], r_config['k']))
    wham_metadata_file.close()
    
if __name__=='__main__':
    main()

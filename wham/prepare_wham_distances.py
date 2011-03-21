#!/usr/bin/python

import os
import sys
import optparse
import numpy
from configobj import ConfigObj
import tempfile
import random

from Queue import Queue, Thread

q = Queue()

def worker():
    while True:
        item = q.get()
        print "Running: "+str(item)
        # run_wham(**item)
        q.task_done()

def run_wham(min=0,max=100,bins=100,metafilepath='',tol=0.0001, temp=315.0):
    pass
    # outfile = "%s_wham.out" % metafilepath
    # command = "wham %f %f %d %f %f 0 %s %s_wham.out" % (min, max, bins, tol, temp, metafilepath, metafilepath)
    # p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # (stdoutdata, stderrdata) = p.communicate()
    # return outfile

def process_config(config_file, start_index=0, end_index=None, output_dir=None, metadata_filename="wham_metadata", percent=100, random=False):
    config = ConfigObj(config_file)
    project_path = os.path.abspath(os.path.dirname(config.filename))

    data_min = None
    data_max = None
    
    if output_dir is None:
        # use a temp dir
        output_dir = tempfile.mkdtemp()
        sys.stderr.write("Creating temporary dir at %s\n" % output_dir)
    elif not output_dir.startswith('/'):
        # relative path
        output_dir = os.path.join(project_path, output_dir)
        if not os.path.exists(output_dir):
            sys.stderr.write("Creating output dir: %s" % output_dir)
            os.mkdir(output_dir)

    # open the wham metadata file for writing
    wham_metadata_file = open(os.path.join(output_dir, metadata_filename), 'w')
    wham_metadata_file.write('# WHAM metadata file generated from Python\n')
    wham_metadata_file.write('# Project Path: %s\n' % project_path)
    wham_metadata_file.write('# Config: %s\n' % config_file)
    wham_metadata_file.write('# Output Dir: %s\n' % output_dir)
    wham_metadata_file.write('# Percent Data Used: %s%%\n' % str(percent))
    wham_metadata_file.write('# Random Selection: %s\n' % str(random))
    
    # loop through replicas in the config file
    for r_id, r_config in config['replicas'].items():
        sys.stderr.write("Extracting data from replica: %s ... " % r_id)
        
        data_file = os.path.join(project_path, r_id, 'distances')
        if not os.path.exists(data_file):
            sys.stderr.write("Data file not found... skipping this replica\n")
            continue
        else:
            field_data = numpy.genfromtxt(data_file)

        data_min = min(field_data.min(), data_min) if data_min else field_data.min()
        data_max = max(field_data.max(), data_max) if data_max else field_data.max()
        
        # first, slice the data
        sample = field_data[start_index:end_index] if end_index else field_data[start_index:] 
        
        # determine the sample size (percent * sample)
        sample_size = int(round(float(percent)/100.0 * float(len(sample))))
        #sys.stderr.write("data/sample: %d/%d\n" % (len(field_data), sample_size))

        if sample_size < 30:
            sys.stderr.write("%s: Sample too small... skipping!\n" % r_id)
            continue
        
        if random:
            final_sample = random.sample(sample, sample_size)
        else:
            final_sample = sample[:sample_size]

        data_file_path = os.path.join(output_dir, 'wham_data_%s' % r_id)
        f = open(data_file_path, 'w')
        for d in final_sample:
            f.write('0.0    %f\n' % d)
        f.close()
        wham_metadata_file.write('%s %s %s\n' % (data_file_path, r_config['coordinate'], str(float(r_config['force']))))
    
    wham_metadata_file.close()
    return {'min':int(data_min), 'max':int(data_max), 'bins': len(config['replicas'].items())*2, 'metafilepath':os.path.join(output_dir, metadata_filename), 'tol':0.0001, 'temp':315.0}

def main():    
    usage = """
        usage: %prog [options] <config.ini> <config.ini> ...

        Prepare and run WHAM.
    """
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-o", "--output-dir", dest="output_dir", default=None, help="Output directory [default: temporary dir]")
    parser.add_option("--combined", dest="combined", default=False, action="store_true", help="Combine all data [default: %default]")
    parser.add_option("--convergence", dest="convergence", default=False, action="store_true", help="Analyze convergence [default: %default]")
    parser.add_option("--error", dest="error", default=False, action="store_true", help="Analyze error [default: %default]")
    parser.add_option("-t", "--threads", dest="worker_threads", default=2, help="Number of WHAM threads to use [default: %default]")
    
    (options, args) = parser.parse_args()
    
    wham_dicts = []
    
    # start the wham threads
    sys.stderr.write("Starting %d worker threads..." % options.worker_threads)
    for i in range(options.worker_threads):
        t = Thread(target=worker)
        t.daemon = True
        t.start()
    
    if options.convergence:
        raise Exception("Convergence not implemented yet")
    elif options.error:
        raise Exception("Error analysis not implemented yet")
    else:
        for config_file in args:
            if not os.path.exists(config_file):
                sys.stderr.write("Config file not found at %s" % config_file)
                continue

            sys.stderr.write("Processing config file: %s" % config_file)
            wham_dict = process_config(config_file)
            wham_dicts.append(wham_dict)
            
            if not options.combined:
                # run it right away
                q.put(wham_dict)
        
        if options.combined:
            # make a single metadatafile containing all the intermediate files
            (fd, fpath) = tempfile.mkstemp()
            sys.stderr.write("Combined metadatafile: %s\n" % fpath)
            outfile = open(fpath, 'w')
            for md in wham_dicts:
                sys.stderr.write("  Combining %s...\n" % md['metafilepath'])
                mfile = open(md['metafilepath'])
                outfile.write(mfile.read())
                mfile.close()
            outfile.close()
            
            combined_dict = wham_dicts[0]
            combined_dict.update({'metafilepath':fpath})
            q.put(combined_dict)
            
    # Wait for wham to finish
    sys.stderr.write("Waiting for WHAM to complete")
    q.join() 

if __name__=='__main__':
    main()

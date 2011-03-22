#!/usr/bin/python

import os
import sys
import optparse
import numpy
from configobj import ConfigObj
import tempfile
import random
import subprocess

from Queue import Queue
from threading import Thread

q = Queue()
outfile_q = Queue()

def combine_metadatas(wham_dicts):
    if len(wham_dicts) == 1:
        return wham_dicts[0]
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
    return combined_dict

def worker():
    while True:
        item = q.get()
        try:
            run_wham(**item)
        except Exception, e:
            sys.stderr.write("\nException while running WHAM: %s\n" % e)
        q.task_done()

def run_wham(min, max, bins, metafilepath, temp=315.0, tol=0.0001, **kwargs):
    outfile = "%s_wham.out" % metafilepath
    sys.stderr.write("Running WHAM: %s -> %s\n" % (metafilepath, outfile))
    command = "wham %f %f %d %f %f 0 %s %s" % (min, max, bins, tol, temp, metafilepath, outfile)
    sys.stderr.write("\t%s\n" % (command))
    p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdoutdata, stderrdata) = p.communicate()
    outfile_q.put(outfile)
    sys.stderr.write("\tDone running WHAM -> %s\n" % (outfile))

def process_pmf(filename, shift=False):
    x = []
    y = []
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
    
    if shift:
        # Calculate the PMF shift by averaging the first 10 bins (bulk water dG)
        s = 1.0*numpy.array(y)[0:10].mean()
        y = [ (i-s) for i in y ]
    
    return zip(x,y)
    
def get_max_bin(pmf):
    maxdG = max(pmf[1])
    return pmf[0][pmf[1].index(maxdG)]

def get_min_bin(pmf):
    mindG = min(pmf[1])
    return pmf[0][pmf[1].index(mindG)]
        
def dG_bind(pmf, imin=0, imax=1):
    """ Calculate the relative binding free energy for a PMF given bounds """
    # dGbind = -(1/kB)
    kB = 0.00198721 # kcal/mol/K
    t = 315.0 # K
    Beta = kB*t
    
    x_prev = None
    dG_prev = None
    s = 0.0
    
    for x,dG in pmf:
        if x > imax:
            break
        elif x > imin and dG_prev is not None:
            # e^{-B*[(dG+dGprev)/2.0]}*(x-x_prev)
            s += numpy.exp((-1.0/Beta)*(dG+dG_prev)/2.0)*(x-x_prev)
        dG_prev = dG
        x_prev = x
    
    dGbind = (-1.0*Beta)*numpy.log(s)
    return dGbind

def process_config(config_file, start_index=0, end_index=None, output_dir=None, metadata_filename="wham_metadata", percent=100, randomize=False):
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
    wham_metadata_file.write('# Random Selection: %s\n' % str(randomize))
    
    # loop through replicas in the config file
    for r_id, r_config in config['replicas'].items():
        sys.stderr.write("%s " % r_id)
        
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
        sys.stderr.write("[%d/%d] " % (len(field_data), sample_size))

        if sample_size < 30:
            sys.stderr.write("\n%s: Sample too small... skipping!\n" % r_id)
            continue
        
        if randomize:
            final_sample = random.sample(field_data, sample_size)
        else:
            final_sample = sample[:sample_size]

        data_file_path = os.path.join(output_dir, 'wham_data_%s' % r_id)
        f = open(data_file_path, 'w')
        for d in final_sample:
            f.write('0.0    %f\n' % d)
        f.close()
        wham_metadata_file.write('%s %s %s\n' % (data_file_path, r_config['coordinate'], str(float(r_config['force']))))
    
    wham_metadata_file.close()
    sys.stderr.write('Done processing config.\n')
    return {'metafilepath':os.path.join(output_dir, metadata_filename)}

def main():    
    usage = """
        usage: %prog [options] <config.ini> <config.ini> ...

        Prepare and run WHAM.
    """
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-o", "--output-dir", dest="output_dir", default=None, help="Output directory [default: temporary dir]")
    parser.add_option("--combined", dest="combined", default=False, action="store_true", help="Combine all data [default: %default]")
    parser.add_option("--convergence", dest="convergence", default=False, action="store_true", help="Analyze convergence [default: %default]")
    parser.add_option("--autoshift", dest="autoshift", default=True, action="store_false", help="Auto-shift PMF [default: %default]")
    parser.add_option("--error", dest="error", default=False, action="store_true", help="Analyze error [default: %default]")
    parser.add_option("--error-blocks", dest="error_blocks", type="int", default=20, help="Number of blocks to sample for error [default: %default]")
    parser.add_option("-t", "--threads", dest="worker_threads", type="int", default=1, help="Number of WHAM threads to use [default: %default]")
    parser.add_option("--wham-min", dest="wham_min", type="float", default=-48, help="Minimum bin value for WHAM [default: %default]")
    parser.add_option("--wham-max", dest="wham_max", type="float", default=0, help="Maximum bin value for WHAM [default: %default]")
    parser.add_option("--wham-bins", dest="wham_bins", type="int", default=200, help="Number of bins for WHAM [default: %default]")
    parser.add_option("--wham-tol", dest="wham_tol", type="float", default=0.0001, help="Tolerance for WHAM [default: %default]")
    parser.add_option("--wham-temp", dest="wham_temp", type="float", default=315.0, help="Temperature for WHAM [default: %default]")
    
    (options, args) = parser.parse_args()
    
    for config_file in args:
        if not os.path.exists(config_file):
            raise Exception("Config file not found at %s\n" % config_file)
    
    # start the wham threads
    sys.stderr.write("Starting %d worker threads...\n" % options.worker_threads)
    for i in range(options.worker_threads):
        t = Thread(target=worker)
        t.daemon = True
        t.start()
            
    wham_defaults = {'min':options.wham_min, 'max':options.wham_max, 'bins':options.wham_bins, 'tol':options.wham_tol, 'temp':options.wham_temp}
    
    if options.convergence:
        # note: ALWAYS COMBINES INPUT CONFIG FILES
        
        # 1) calculate blocks of data in sequential order for each config file
        #       block size is 10% of the max n
        # 2) calculate the PMFs from each block
        # 3) calculate some value from each PMF (dG_bind)
        # 4) print <block>,<value> to plot
        
        
        # store the max N for each config file
        max_n = {}
        # store the current block index for each config file
        current_block = {}
        
        
        raise Exception("Convergence not implemented yet")
    elif options.error:
        # note: ALWAYS COMBINES INPUT CONFIG FILES
        
        # 1) calculate random blocks of data for each config file
        # 2) calculate PMFs for each block
        # 3) plot each PMF, get max/min values per bin, stdev per bin
        
        i=0        
        while i < options.error_blocks:
            wham_dicts = []
            for config_file in args:
                sys.stderr.write("Processing config file: %s\n" % config_file)                
                md = process_config(config_file, percent=25, randomize=True)
                md.update(wham_defaults)
                wham_dicts.append(md)
            combined_dict = combine_metadatas(wham_dicts)
            q.put(combined_dict)
            i += 1
                
        # Wait for wham to finish
        sys.stderr.write("Waiting for WHAM to complete\n")
        q.join()
        
        # process the PMFs
        error_data = {}
        outfile = outfile_q.get_nowait()
        while True:
            sys.stderr.write("Processing WHAM outfile: %s\n" % outfile)
            pmf = process_pmf(outfile, shift=options.autoshift)
            
            print dG_bind(pmf, imin=-31.0, imax=-21.0)
            
            for x,y in pmf:
                if x not in error_data:
                    error_data[x] = [y]
                else:
                    error_data[x].append(y)
            
            outfile_q.task_done()
            try:
                outfile = outfile_q.get_nowait()
            except:
                sys.stderr.write("All outfiles processed...\n")
                break

        # now we have the combined PMF data
        (fd, fpath) = tempfile.mkstemp()
        outfile = open(fpath, 'w')
        sys.stderr.write("Writing errors\n")
        outfile.write('BIN,MEAN,SEM,STDEV,MIN,MAX\n')
        
        for key in sorted(error_data.keys()):
            d = numpy.array(error_data[key])
            sys.stderr.write("%0.2f,%d " % (key, d.size))
            sem = numpy.std(d, ddof=1)/numpy.sqrt(d.size)
            outfile.write('%f,%f,%f,%f,%f,%f\n' % (key,d.mean(),sem,numpy.std(d),d.min(),d.max()))
        outfile.close()
        sys.stderr.write("\n"+fpath+"\n")
        
    else:
        # standard procedure
        # TODO: remove combined option?
        wham_dicts = []
        for config_file in args:
            sys.stderr.write("Processing config file: %s\n" % config_file)
            wham_dict = process_config(config_file)
            wham_dict.update(wham_defaults)
            wham_dicts.append(wham_dict)
            
            if not options.combined:
                # run it right away
                q.put(wham_dict)
        
        if options.combined:
            combined_dict = combine_metadatas(wham_dicts)
            q.put(combined_dict)
            
        # Wait for wham to finish
        sys.stderr.write("Waiting for WHAM to complete\n")
        q.join() 
    
        try:
            outfile = outfile_q.get_nowait()
            while True:
                print outfile
                outfile_q.task_done()
                outfile = outfile_q.get_nowait()
        except:
            # queue empty
            pass


if __name__=='__main__':
    main()

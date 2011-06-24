#!/usr/bin/python

import os
import sys
import optparse
from configobj import ConfigObj
import shutil
import subprocess

from Queue import Queue
from threading import Thread

q = Queue()

def gunzip(gzfile):
    command = "/usr/bin/gunzip %s" % (gzfile)
    p = subprocess.Popen(command, shell=True)
    p.wait()
    return gzfile.replace('.gz','')

def VMD_filter_pdb(pdb_file, psf_file, selection="protein", outfile="outfile.pdb"):
    tcl = """set A [atomselect top "%s"]
$A writepdb %s
quit
    """ % (selection, outfile)
    command = "vmd -dispdev none -psf %s -pdb %s" % (psf_file, pdb_file)
    p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdoutdata, stderrdata) = p.communicate(input=tcl)

def worker():
    while True:
        item = q.get()
        try:
            sys.stderr.write("Processing: %s\n" % item['pdb'])
            # run VMD on the pdb file
            pdb = gunzip(item['pdb'])
            VMD_filter_pdb(pdb, item['options'].psf_file, selection="protein or ions", output_filename=pdb+'_')
            shutil.move(pdb+'_',pdb)
        except Exception, e:
            sys.stderr.write("\nException while running VMD: %s\n" % e)
        q.task_done()

def main():    
    usage = """
        usage: %prog [options] <config.ini>
    """
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-o", "--output-dir", dest="output_dir", default='/dev/shm', help="Output directory [default: %default]")
    parser.add_option("-t", "--threads", dest="worker_threads", type="int", default=1, help="Number of WHAM threads to use [default: %default]")
    # parser.add_option("-s", "--selection", dest="selection", default="protein", help="VMD selection to keep [default: %default]")
    parser.add_option("-p", "--psf", dest="psf_file", default=None, help="PSF structure file [default: %default]")
    
    (options, args) = parser.parse_args()
    
    config_file = args[0]
    
    if not os.path.exists(config_file):
        raise Exception("Config file not found at %s\n" % config_file)
    
    # start the wham threads
    sys.stderr.write("Starting %d worker threads...\n" % options.worker_threads)
    for i in range(options.worker_threads):
        t = Thread(target=worker)
        t.daemon = True
        t.start()
    
    config = ConfigObj(config_file)
    project_path = os.path.abspath(os.path.dirname(config.filename))
    scratch_path = project_path.replace('/project/pomes/dacaplan/','/scratch/dacaplan/scratch/')
    sys.stderr.write("Scratch path: %s\n" % scratch_path)

    reps = []

    # get the replicas from the config by ID and COORD
    for r_id, r_config in config['replicas'].items():
        c = float(r_config['coordinate'])
        reps.append((r_id, c))

    # for each replica (sorted by coord) copy all the PDBs to the output dir
    i=0
    for r in sorted(reps, key=lambda r: r[1]):
        pdb_path = os.path.join(scratch_path, r[0])
        for dirname, dirnames, filenames in os.walk(pdb_path):
            for filename in filenames:
                if filename.endswith('.pdb.gz'):
                    pdb = os.path.join(options.output_dir, '%d.pdb.gz' % i)
                    shutil.copy(os.path.join(dirname, filename), pdb)
                    # send it to be processed by VMD
                    q.put({'options': options, 'pdb':pdb})
                    i = i+1
                    break
            break
        break
                        
    sys.stderr.write("Waiting for jobs to complete\n")
    q.join()
    
    # now run catdcd on all the PDBs
    
    
    # try:
    #     outfile = outfile_q.get_nowait()
    #     while True:
    #         print outfile
    #         outfile_q.task_done()
    #         outfile = outfile_q.get_nowait()
    # except:
    #     # queue empty
    #     pass


if __name__=='__main__':
    main()

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

SEL='segname PEPA or ions'

def gunzip(gzfile):
    command = "/usr/bin/gunzip -f %s" % (gzfile)
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
    return p.communicate(input=tcl)


def VMD_write_psf(pdb_file, psf_file, selection="protein", outfile="outfile.psf"):
    tcl = """set A [atomselect top "%s"]
$A writepsf %s
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
            pdb = gunzip(item['pdb'])

            # write new PSF
            if item['write_psf']:
                sys.stderr.write("Writing PSF file\n")
                newpsf = item['options'].psf_file.replace('.psf','_trim.psf')
                VMD_write_psf(pdb, item['options'].psf_file, selection=SEL, outfile=newpsf)
                sys.stderr.write("Wrote trimmed PSF file: %s\n" % newpsf)

            # run VMD on the pdb file
            (stdout, stderr) = VMD_filter_pdb(pdb, item['options'].psf_file, selection=SEL, outfile=pdb+'_')
            print stderr
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
    parser.add_option("-f", "--frames", dest="frames_per_replica", type="int", default=-1, help="Frames per replica [default: all frames]")
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
    write_psf = True
    for r in sorted(reps, key=lambda r: r[1]):
        pdb_path = os.path.join(scratch_path, r[0])
        for dirname, dirnames, filenames in os.walk(pdb_path):
            for j,filename in enumerate(filenames):
                if j == options.frames_per_replica:
                    break
                if filename.endswith('.pdb.gz'):
                    if not os.path.exists(os.path.join(options.output_dir,'%d.pdb'%i)):
                        pdb = os.path.join(options.output_dir, '%d.pdb.gz' % i)
                        shutil.copy(os.path.join(dirname, filename), pdb)
                        # send it to be processed by VMD
                        # HACK: this is stupid
                        if write_psf:
                            q.put({'options': options, 'pdb':pdb, 'write_psf': True})
                            write_psf = False
                        else:
                            q.put({'options': options, 'pdb':pdb, 'write_psf': False})
                    else:
                        sys.stderr.write("Skipping PDB %d\n" % i)
                    i = i+1
                        
    sys.stderr.write("Waiting for jobs to complete\n")
    q.join()
    
    sys.stderr.write("Running catdcd\n")
    # now run catdcd on all the PDBs
    #newpsf = options.psf_file.replace('.psf','_trim.psf')
    #cmd = "catdcd -o /dev/shm/%s.dcd -stype psf -s %s %s" % (config['title'], newpsf, ' '.join([ '-pdb %s' % os.path.join(options.output_dir, '%d.pdb' % j) for j in range(i) ]))
    #print "catdcd -o out.dcd -s %s %s" % (newpsf, ' '.join([ '-pdb %d.pdb' % j for j in range(i) ]))
    
    cmd = "cd %s;catdcd -o %s.dcd -s 0.pdb `ls *.pdb | sort -n | xargs -s 50 -I %% echo -n \"-pdb %% \"`" % (options.output_dir, config['title'])
    print cmd
    
    p = subprocess.Popen(cmd, shell=True)
    p.wait()
    
if __name__=='__main__':
    main()

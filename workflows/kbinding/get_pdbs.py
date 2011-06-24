#!/usr/bin/python

import os
import shutil
import sys
import optparse
import numpy
from configobj import ConfigObj

def main():    
    usage = """
        usage: %prog [options] <config.ini>
    """
    
    parser = optparse.OptionParser(usage)
    (options, args) = parser.parse_args()
    
    config = ConfigObj(args[0])
    project_path = os.path.abspath(os.path.dirname(config.filename))
    scratch_path = project_path.replace('/project/pomes/dacaplan/','/scratch/dacaplan/scratch/')
    sys.stderr.write("Scratch path: %s\n" % scratch_path)

    reps = []

    for r_id, r_config in config['replicas'].items():
        c = float(r_config['coordinate'])
        reps.append((r_id, c))
    
    i=0
    for r in sorted(reps, key=lambda r: r[1]):
        # get list of *pdb* 
        pdb_path = os.path.join(scratch_path, r[0])
        for dirname, dirnames, filenames in os.walk(pdb_path):
            for filename in filenames:
                if filename.endswith('pdb'):
                    shutil.copy(os.path.join(dirname, filename), '/dev/shm/%d.pdb' % i)
                    i = i+1
                elif filename.endswith('gz'):
                    shutil.copy(os.path.join(dirname, filename), '/dev/shm/%d.pdb.gz' % i)
                    i = i+1
        
if __name__=='__main__':
    main()

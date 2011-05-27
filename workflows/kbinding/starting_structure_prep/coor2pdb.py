#!/usr/bin/python

import os
import optparse
import subprocess

def main():    
    usage = """
        usage: %prog [options] <PSF> <COOR>
    """
    
    parser = optparse.OptionParser(usage)
    (options, args) = parser.parse_args()
    
    if len(args) == 0:
        parser.error("Please specify an input file")
    elif not os.path.exists(args[0]):
        parser.error("Input file %s not found!" % args[0])

    psf = args[0]
    coor = args[1]
    
    tcl = """set A [atomselect top "all"]
$A writepdb coor_md.pdb
quit
    """

    command = "vmd -dispdev none -psf %s -namdbin %s" % (psf, coor)
    p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdoutdata, stderrdata) = p.communicate(input=tcl)

if __name__=='__main__':
    main()




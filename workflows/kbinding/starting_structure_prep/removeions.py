#!/usr/bin/python

import os
import optparse
import subprocess

VMD_PATH = "vmd"

def main():    
    usage = """
        usage: %prog [options] <PSF> <PDB/DCD>
    """
    
    parser = optparse.OptionParser(usage)
    (options, args) = parser.parse_args()
    
    if len(args) == 0:
        parser.error("Please specify an input file")
    elif not os.path.exists(args[0]):
        parser.error("Input file %s not found!" % args[0])

    psf = args[0]
    pdb = args[1]
    
    tcl = """set A [atomselect top "all and not (name POT and resid > 4719) and not (name CLA)"]
$A writepsf md.psf
$A writepdb md.pdb
quit
    """

    command = "vmd -dispdev none -psf %s -pdb %s" % (psf, pdb)
    p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdoutdata, stderrdata) = p.communicate(input=tcl)

if __name__=='__main__':
    main()




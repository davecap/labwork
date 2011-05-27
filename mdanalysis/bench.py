#!/usr/bin/python
# -*- coding: utf-8 -*-

import optparse
from MDAnalysis import * # analysis of NAMD trajectories

def main():
    usage = """
        usage: %prog [options] <PSF file> <PDB file> <DCD file>
    """
    
    # initialize the parser with our custom usage string (above)
    parser = optparse.OptionParser(usage)
    # parse the command line
    (options, args) = parser.parse_args()
    
    # check to make sure all arguments are there
    if len(args) < 2:
        parser.error("No input files specified")
    
    psf_file = args[0]
    pdb_file = args[1]
    
    # load the PSF and PDB into the analysis system
    print "Loading reference system: %s, %s" % (psf_file, pdb_file)
    ref = Universe(psf_file, pdb_file, permissive=True)

    try:
        dcd_file = args[2]
    except:
        pass
    else:
        print "Loading trajectory: %s" % (dcd_file)
        trj = Universe(psf_file, dcd_file)

if __name__ == '__main__':
    main()

#!/usr/bin/python
# -*- coding: utf-8 -*-

# import required modules
import optparse # command-line parsing
import numpy # numerical operations
from tables import * # pytables (database)
from MDAnalysis import * # analysis of NAMD trajectories

# my analysis modules (in the same directory as this file)
from analysis import *
from rmsd import RMSD
from hbonds import HydrogenBondAnalysis
from nearby import NearbyCountAnalysis

# main() isn't special, it is just like any other Python function
def main():
    usage = """
        usage: %prog [options] <PSF file> <PDB file> <DCD file> <DCD file 1> <DCD file 2>
    """
    
    # initialize the parser with our custom usage string (above)
    parser = optparse.OptionParser(usage)
    # parse the command line
    (options, args) = parser.parse_args()
    
    # check to make sure all arguments are there
    if len(args) < 3:
        parser.error("No input files specified")
    
    psf_file = args[0]
    pdb_file = args[1]
    
    # load the PSF and PDB into the analysis system
    print "Loading reference system: %s, %s" % (psf_file, pdb_file)
    ref = Universe(psf_file, pdb_file)
        
    n = NearbyCountAnalysis('segid PEPA and resid 139', 'resname POT', cutoff=3.5)
    n1 = NearbyCountAnalysis('segid PEPA and resid 132', 'resname POT', cutoff=3.5)

    for dcd_file in args[2:]:
        print "Loading trajectory: %s" % (dcd_file)
        trj = Universe(psf_file, dcd_file)
        n.prepare(trj=trj)
        n1.prepare(trj=trj)

        for ts in trj.trajectory:
            res = n.process(ts.frame)
            if res:
                print "POT within 3.5 A of 139:"
                print res
                print ts
            res = n1.process(ts.frame)
            if res:
                print "POT within 3.5 A of 132:"
                print res
                print ts
        

if __name__ == '__main__':
    main()

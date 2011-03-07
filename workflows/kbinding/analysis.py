#!/usr/bin/python

import os
import optparse
#from configobj import ConfigObj, flatten_errors
import subprocess

VMD_PATH = "vmd"

def main():    
    usage = """
        usage: %prog [options] <PDB/DCD file>
    """
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-p", "--psf", dest="psf_file", default=None, help="PSF structure file for DCDs [default: %default]")
    
    (options, args) = parser.parse_args()
    
    if len(args) == 0:
        parser.error("Please specify an input file")
    elif not os.path.exists(args[0]):
        parser.error("Input file %s not found!" % args[0])
    
    if args[0].lower().endswith('.dcd'):
        if options.psf_file is None:
            parser.error("PSF file required (-p) if using a DCD trajectory file")
        VMD_analyze(psf_file=options.psf_file, dcd_file=args[0])   
    elif args[0].lower().endswith('.pdb'):
        VMD_analyze(psf_file=options.psf_file, pdb_file=args[0])   
    else:
        parser.error("Invalid input file, must be a PDB or DCD file")

def VMD_analyze(psf_file, dcd_file=None, pdb_file=None):
    tcl = """
set e286ca [atomselect top "segname PEPA and resid 286 and name CA"]
set k4712 [atomselect top "segname POT and resid 4712"]

set num_steps [molinfo 0 get numframes]
for {set frame 0} {$frame < $num_steps} {incr frame} {
    set e286caz [$e286ca get {z}]
    set k4712z [$k4712 get {z}]
    puts "PARSE [ expr $k4712z-$e286caz ]"
}
quit
    """

    if dcd_file:
        command = "vmd -dispdev none -psf %s -dcd %s" % (psf_file, dcd_file)
    elif pdb_file:
        command = "vmd -dispdev none -psf %s -pdb %s" % (psf_file, pdb_file)
    p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdoutdata, stderrdata) = p.communicate(input=tcl)
    
    #print stdoutdata
    #print stderrdata
    for line in stdoutdata.split('\n'):
        if line.startswith('PARSE'):
            line_parts = line.split(' ')
            for d in line_parts[1:]:
                print d
 
if __name__=='__main__':
    main()




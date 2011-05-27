#!/usr/bin/python

import os
import optparse
import sys
#from configobj import ConfigObj, flatten_errors
import subprocess

VMD_PATH = "vmd"

def main():    
    usage = """
        usage: %prog [options] <PDB/DCD file>
    """
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-s", "--pdb", dest="pdb_file", default=None, help="PDB file [default: %default]")
    parser.add_option("-p", "--psf", dest="psf_file", default=None, help="PSF structure file for DCDs [default: %default]")
    
    (options, args) = parser.parse_args()
    
    if len(args) == 0:
        parser.error("Please specify an input file")
    elif not os.path.exists(args[0]):
        parser.error("Input file %s not found!" % args[0])
    
    if args[0].lower().endswith('.dcd'):
        if options.psf_file is None:
            parser.error("PSF file required (-p) if using a DCD trajectory file")
        VMD_analyze_dcd(pdb_file=options.pdb_file, psf_file=options.psf_file, dcd_file=args[0])   

def VMD_analyze_dcd(pdb_file, psf_file, dcd_file):
    tcl = """
set mol [mol new %s type psf waitfor all]
set all [atomselect $mol all]
set ref [atomselect $mol "name CA" frame 0]
set sel [atomselect $mol "name CA"]
set wat [atomselect $mol "(resname TIP3 and name OH2) and within 5 of protein"]
mol addfile %s type pdb waitfor all
mol addfile %s type dcd waitfor all

set num_steps [molinfo 0 get numframes]
for {set frame 0} {$frame < $num_steps} {incr frame} {
    $all frame $frame
    $sel frame $frame
    $ref frame $frame
    $all move [measure fit $sel $ref]
}

volmap density $wat -res 0.5 -weight mass -allframes -combine avg -mol top -checkpoint 0 -o %s_volmap.dx
quit
    """ % (psf_file, pdb_file, dcd_file, dcd_file)
    sys.stderr.write(tcl+'\n')

    command = "vmd -dispdev none"
    p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdoutdata, stderrdata) = p.communicate(input=tcl)

    #sys.stderr.write(stderrdata+'\n')
    print stdoutdata
    print stderrdata

    # for line in stdoutdata.split('\n'):                                                                                                                                                                                            
    #     if line.startswith('RMSD'):                                                                                                                                                                                           
    #         line_parts = line.split(' ')                                                                                                                                                                                           
    #         for d in line_parts[1:]:                                                                                                                                                                                               
    #             print d    
    
if __name__=='__main__':
    main()




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
    parser.add_option("-o", "--output-prefix", dest="output_prefix", default="filtered", help="Output file prefix [default: %default]")
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
        VMD_analyze_dcd(pdb_file=options.pdb_file, psf_file=options.psf_file, dcd_files=args, output_prefix=options.output_prefix)   

def VMD_analyze_dcd(pdb_file, psf_file, dcd_files, output_prefix="filtered"):
    output_psf = "%s.psf" % output_prefix
    output_dcd = "%s.dcd" % output_prefix
        
    tcl = """
source bigdcd.tcl
proc myrmsd { frame } {
    global ref sel all
    $all move [measure fit $sel $ref]
    puts "RMSD $frame,[measure rmsd $sel $ref]"
}
set mol [mol new %s type psf waitfor all]
set all [atomselect $mol all]
set ref [atomselect $mol "name CA" frame 0]
set sel [atomselect $mol "name CA"]
mol addfile %s type pdb waitfor all
bigdcd myrmsd %s
bigdcd_wait
quit
    """ % (psf_file, pdb_file,' '.join(dcd_files))
    sys.stderr.write(tcl+'\n')

    command = "vmd -dispdev none"
    p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdoutdata, stderrdata) = p.communicate(input=tcl)

    sys.stderr.write(stderrdata+'\n')
    #print stdoutdata
    #print stderrdata

    for line in stdoutdata.split('\n'):                                                                                                                                                                                            
        if line.startswith('RMSD'):                                                                                                                                                                                           
            line_parts = line.split(' ')                                                                                                                                                                                           
            for d in line_parts[1:]:                                                                                                                                                                                               
                print d    
    
if __name__=='__main__':
    main()




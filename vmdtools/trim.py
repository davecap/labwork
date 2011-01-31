#!/usr/bin/python

import os
import optparse
from configobj import ConfigObj, flatten_errors
import subprocess

VMD_PATH = "vmd"

def main():    
    usage = """
        usage: %prog [options] <PDB/DCD file>
    """
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-o", "--output-prefix", dest="output_prefix", default="filtered", help="Output file prefix [default: %default]")
    parser.add_option("-s", "--selection", dest="selection", default="protein", help="VMD selection to keep [default: %default]")
    parser.add_option("-p", "--psf", dest="psf_file", default=None, help="PSF structure file for DCDs [default: %default]")
    
    (options, args) = parser.parse_args()
    
    if len(args) == 0:
        parser.error("Please specify an input file")
    elif not os.path.exists(args[0]):
        parser.error("Input file %s not found!" % args[0])
    
    if args[0].lower().endswith('.dcd'):
        if options.psf_file is None:
            parser.error("PSF file required (-p) if using a DCD trajectory file")
        VMD_filter_dcd(psf_file=options.psf_file, dcd_file=args[0], selection=options.selection, output_prefix=options.output_prefix)   
    elif args[0].lower().endswith('.pdb'):
        VMD_filter_pdb(args[0], selection=options.selection, output_prefix=options.output_prefix)
    else:
        parser.error("Invalid input file, must be a PDB or DCD file")

def VMD_filter_dcd(psf_file, dcd_file, selection="protein", output_prefix="filtered"):
    output_psf = "%s.psf" % output_prefix
    output_dcd = "%s.dcd" % output_prefix
        
    tcl = """set A [atomselect top "%s"]
    $A writepsf %s
    animate write dcd %s sel $A    
    quit
    """ % (selection, output_psf, output_dcd)

    command = "vmd -dispdev none -psf %s -dcd %s" % (psf_file, dcd_file)
    p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdoutdata, stderrdata) = p.communicate(input=tcl)
 
def VMD_filter_pdb(pdb_file, selection="protein", output_prefix="filtered"):
    output_pdb = "%s.pdb" % output_prefix
    
    tcl_template = """set A [atomselect top "%s"]
    $A writepdb %s
    quit
    """ % (selection, output_pdb)
    
    command = "vmd -dispdev none -pdb %s" % (pdb_file)   
    p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdoutdata, stderrdata) = p.communicate(input=tcl)

if __name__=='__main__':
    main()




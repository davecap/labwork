#!/usr/bin/python

import os
import optparse
from configobj import ConfigObj, flatten_errors
import subprocess

VMD_PATH = "vmd"

def main():    
    usage = """
        Align a PDB of DCD trajectory to a structure template
        usage: %prog [options] <PDB/DCD file> <PDB file template>
    """
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-o", "--output-prefix", dest="output_prefix", default="aligned", help="Output file prefix [default: %default]")
    parser.add_option("-s", "--selection", dest="selection", default="protein and backbone", help="VMD selection to align by3 [default: %default]")
    parser.add_option("-p", "--psf", dest="psf_file", default=None, help="PSF structure file for DCDs [default: %default]")
    
    (options, args) = parser.parse_args()
    
    if len(args) < 2:
        parser.error("Please specify the structure and template files")
    elif not os.path.exists(args[0]):
        parser.error("Input file %s not found!" % args[0])
    elif not os.path.exists(args[1]):
        parser.error("Input file %s not found!" % args[1])
    
    if args[0].lower().endswith('.dcd'):
        if options.psf_file is None:
            parser.error("PSF file required (-p) if using a DCD trajectory file")
        VMD_align_dcd(psf_file=options.psf_file, template_pdb=args[1], selection=options.selection, dcd_file=args[0], output_prefix=options.output_prefix)   
    elif args[0].lower().endswith('.pdb'):
        VMD_align_pdb(args[0], template_pdb=args[1], selection=options.selection, output_prefix=options.output_prefix)
    else:
        parser.error("Invalid input file, must be a PDB or DCD file")

def VMD_align_dcd(psf_file, dcd_file, template_pdb, selection, output_prefix="filtered"):
    output_dcd = "%s.dcd" % output_prefix
    
    tcl = """mol load pdb "%s"
set sel0all [atomselect 0 all]
set sel0 [atomselect 0 "%s"]
set sel1 [atomselect 1 "%s"]

set num_steps [molinfo 0 get numframes]
for {set frame 0} {$frame < $num_steps} {incr frame} {
    $sel0 frame $frame
    $sel0all frame $frame
    #set rmsd [measure rmsd $sel0 $sel1]
    #puts "RMSD before alignment: $rmsd"
    set trans_mat [measure fit $sel0 $sel1]
    $sel0all move $trans_mat
    #set rmsd [measure rmsd $sel0 $sel1]
    #puts "RMSD after alignment: $rmsd"
}
mol top 0
set sel0all [atomselect 0 all]
animate write dcd "%s" sel $sel0all
quit
    """ % (template_pdb, selection, selection, output_dcd)

    command = "vmd -dispdev none -psf %s -dcd %s" % (psf_file, dcd_file)
    p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdoutdata, stderrdata) = p.communicate(input=tcl)
 
def VMD_align_pdb(pdb_file, template_pdb, selection, output_prefix="filtered"):
    output_pdb = "%s.pdb" % output_prefix
    
    tcl = """mol load pdb "%s"
set sel0all [atomselect 0 all]
set sel0 [atomselect 0 "%s"]
set sel1 [atomselect 1 "%s"]
#set rmsd [measure rmsd $sel0 $sel1]
#puts "RMSD before alignment: $rmsd"
set trans_mat [measure fit $sel0 $sel1]
$sel0 move $trans_mat
#set rmsd [measure rmsd $sel0 $sel1]
#puts "RMSD after alignment: $rmsd"
$sel0all writepdb "%s"
quit
    """ % (template_pdb, selection, selection, output_pdb)
    
    command = "vmd -dispdev none -pdb %s" % (pdb_file)
    p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdoutdata, stderrdata) = p.communicate(input=tcl)

if __name__=='__main__':
    main()




#!/usr/bin/python
"""
Stage: Prepare a starting structure.

Stage takes 2 or 3 input files:
    1) Staging structure (PDB or DCD)
    2) Starting structure (PDB)
    3) (optional) PSF file if DCD is used
    
The user specifies a selection string, which represents the atoms to be staged.
First, an alignment of each of the staging structures is carried out against the starting structure.
Then, the selected atoms in the starting structure are moved to their locations in the staging structure(s).

One PDB is output for each staging structure provided.
"""

import os
import optparse
from configobj import ConfigObj, flatten_errors
import subprocess
from string import Template

VMD_PATH = "vmd"

def main():    
    usage = """
        Stage: Prepare a starting structure.
        
        usage: %prog [options] <staging structure/traj> <starting structure>
    """
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-o", "--output-prefix", dest="output_prefix", default="aligned", help="Output file prefix [default: %default]")
    parser.add_option("-a", "--align-selection", dest="align_selection", default="protein and backbone", help="VMD selection to align by [default: %default]")
    parser.add_option("-t", "--target-selection", dest="target_selection", default="resid 1", help="VMD target selection to stage [default: %default]")
    parser.add_option("-p", "--psf", dest="psf_file", default=None, help="PSF structure file for DCDs [default: %default]")
    parser.add_option("-v", "--verbose", dest="verbose", default=False, action="store_true", help="Verbose mode! [default: %default]")
    
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
        tcl = dcd_stage_tcl(staging_psf=options.psf_file, staging_dcd=args[0], starting_pdb=args[1], align_selection=options.align_selection, target_selection=options.target_selection, output_prefix=options.output_prefix)   
    elif args[0].lower().endswith('.pdb'):
        tcl = pdb_stage_tcl(args[0], starting_pdb=args[1], align_selection=options.align_selection, target_selection=options.target_selection, output_prefix=options.output_prefix)
    else:
        parser.error("Invalid input file, must be a PDB or DCD file")
    
    run_VMD(tcl, options.verbose)

def dcd_stage_tcl(staging_psf, staging_dcd, starting_pdb, align_selection, target_selection, output_prefix="staged"):
    data = {    'staging_psf': staging_psf,
                'staging_dcd': staging_dcd,
                'starting_pdb': starting_pdb,
                'align_selection': align_selection,
                'target_selection': target_selection,
                'output_prefix': output_prefix,
            }

    tcl_template = """mol load psf "$staging_psf" dcd "$staging_dcd"
mol load pdb "$starting_pdb"
# staging: 0
# starting: 1

# align staging -> starting (0->1)
set sel0all [atomselect 0 all]
set sel1all [atomselect 1 all]
set sel0align [atomselect 0 "$align_selection"]
set sel1align [atomselect 1 "$align_selection"]
set sel0target [atomselect 0 "$target_selection"]
set sel1target [atomselect 1 "$target_selection"]

set num_steps [molinfo 0 get numframes]
for {set frame 0} {$$frame < $$num_steps} {incr frame} {
    $$sel0align frame $$frame
    $$sel0all frame $$frame
    $$sel0target frame $$frame
    set trans_mat [measure fit $$sel0align $$sel1align]
    $$sel0all move $$trans_mat
    # move target in starting (1) to position of the same target in staging (0)
    set trans_mat [measure fit $$sel1target $$sel0target]
    $$sel1target move $$trans_mat
    $$sel1all writepdb "${output_prefix}_$$frame.pdb"
}
quit
"""
    return Template(tcl_template).substitute(data)
 
def pdb_stage_tcl(staging_pdb, starting_pdb, align_selection, target_selection, output_prefix="staged"):    
    data = {    'staging_pdb': staging_pdb,
                'starting_pdb': starting_pdb,
                'align_selection': align_selection,
                'target_selection': target_selection,
                'output_prefix': output_prefix,
            }
    
    tcl_template = """mol load pdb "$staging_pdb"
mol load pdb "$starting_pdb"
# staging: 0
# starting: 1

# align staging -> starting (0->1)
set sel0all [atomselect 0 all]
set sel1all [atomselect 1 all]
set sel0align [atomselect 0 "$align_selection"]
set sel1align [atomselect 1 "$align_selection"]
set trans_mat [measure fit $$sel0align $$sel1align]
$$sel0all move $$trans_mat

# move target in starting (1) to position of the same target in staging (0)
set sel0target [atomselect 0 "$target_selection"]
set sel1target [atomselect 1 "$target_selection"]
set trans_mat [measure fit $$sel1target $$sel0target]
$$sel1target move $$trans_mat

$$sel1all writepdb "$output_prefix.pdb"
quit
""" 
    return Template(tcl_template).substitute(data)    

def run_VMD(tcl, verbose=False):
    command = "%s -dispdev none" % (VMD_PATH)
    
    if verbose:
        print command
        print tcl
    p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdoutdata, stderrdata) = p.communicate(input=tcl)
    if verbose:
        print stdoutdata
        print stderrdata

if __name__=='__main__':
    main()




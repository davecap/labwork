import optparse
import subprocess
from string import Template

def pdb_stage_tcl(staging_pdb, starting_pdb, align_selection, template_target_selection, starting_target_selection, output_pdb="staged"):    
    data = {    'staging_pdb': staging_pdb,
                'starting_pdb': starting_pdb,
                'align_selection': align_selection,
                'template_target_selection': template_target_selection,
                'starting_target_selection': starting_target_selection,
                'output_pdb': output_pdb,
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
set sel0target [atomselect 0 "$template_target_selection"]
set sel1target [atomselect 1 "$starting_target_selection"]
set trans_mat [measure fit $$sel1target $$sel0target]
$$sel1target move $$trans_mat

$$sel1all writepdb "$output_pdb"
quit
""" 
    return Template(tcl_template).substitute(data)    

def run_VMD(tcl):
    command = "vmd -dispdev none"
    p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return p.communicate(input=tcl)

def main():    
    usage = """
        usage: %prog [options] <starting structure PDB> <template PDB> <output PDB>
    """
    parser = optparse.OptionParser(usage)
    parser.add_option("-t", dest="template_target", default='segname POT and resid 4712', help="Template target selection")
    parser.add_option("-s", dest="starting_target", default='segname POT and resid 4719', help="Starting target selection")

    options, args = parser.parse_args()
    
    try:
        starting_pdb = args[0]
        template_pdb = args[1]
        output_pdb = args[2]
    except:
        parser.error("Please verify parameters")
    
    # Align using the VMD TCL script
    tcl = pdb_stage_tcl(staging_pdb=template_pdb, starting_pdb=starting_pdb, align_selection="protein and backbone", template_target_selection=options.template_target, starting_target_selection=options.starting_target, output_pdb=output_pdb)
    (stderr, stdout) = run_VMD(tcl)
    #print stderr
    #print stdout
                    
if __name__ == "__main__":
    main()






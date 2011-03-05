import optparse
import sys
import os

import umbrellas
import stage

def main():    
    usage = """
        usage: %prog [options] <template_config.ini> <config.ini> <starting structure PDB>
    """
    parser = optparse.OptionParser(usage)
    options, args = parser.parse_args()
    
    try:
        template_config = args[0]
        output_config = args[1]
        starting_pdb = args[2]
    except:
        parser.error("Please verify parameters")
    
    template_ensemble = umbrellas.Ensemble(template_config)
    output_ensemble = umbrellas.Ensemble(output_config)
    
    output_path = os.path.dirname(output_config)
    output_pdb_path = os.path.dirname(output_config)+'/pdbs'
    os.mkdir(output_pdb_path)
        
    # for each template replica
    #   align a with starting structure
    #   create new output replica and PDB
    for i,t in enumerate(template_ensemble.get_replicas()):
        template_pdb = t.get_parameter('coordinates')
        template_psf = template_ensemble.config['psf']
        output_pdb = output_pdb_path+'/%d.pdb'%i
        
        tcl = stage.pdb_stage_tcl(staging_pdb=template_pdb, starting_pdb=starting_pdb, align_selection="protein and backbone", target_selection="protein and backbone", output_pdb=output_pdb)
        #run_VMD(tcl, options.verbose)
        print tcl
        
        #r = output_ensemble.add_replica(name='%d'%i, **t._parameters)
        #r._parameters['coordinate'] = output_pdb
                    
if __name__ == "__main__":
    main()

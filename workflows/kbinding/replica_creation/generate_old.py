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
    
    template_path = template_ensemble.basedir()
    output_path = output_ensemble.basedir()
    output_pdb_path = os.path.join(output_path,'pdbs')

    if not os.path.exists(output_pdb_path):
        sys.stderr.write("Creating output PDB path: %s\n" % output_pdb_path)
        os.mkdir(output_pdb_path)
        
    # for each template replica
    #   align a with starting structure
    #   create new output replica and PDB
    for i,t in enumerate(template_ensemble.get_replicas()):
        sys.stderr.write("Aligning template %d\n" % i)
        template_pdb = os.path.join(template_path,t.parameter('coordinates'))
        template_psf = os.path.join(output_path, template_ensemble.config['psf'])
        output_pdb = os.path.join(output_pdb_path,'%d.pdb'%i)
        
        # Align using the VMD TCL script
        tcl = stage.pdb_stage_tcl(staging_pdb=template_pdb, starting_pdb=starting_pdb, align_selection="protein and backbone", target_selection="protein and backbone", output_pdb=output_pdb)
        stage.run_VMD(tcl)
        
        # Create the new replica with the same parameters but set the PDB path
        r = output_ensemble.add_replica(name='%d'%i, **t.parameters)
        r.parameters['coordinates'] = 'pdbs/%d.pdb'%i
        output_ensemble.save()

    sys.stderr.write("Created %d replicas from starting structure %s" % (i+1, starting_pdb))
                    
if __name__ == "__main__":
    main()

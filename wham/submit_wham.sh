#!/bin/sh

cd /hpf/projects1/pomes/dacaplan/projects/cytoxdrep/wildtype/all_subunits/chi_one_pmf/pmf_e_cancel/dihedrals
/hpf/projects1/pomes/rhenry/exe/wham-release-2.0.1/wham-2d/wham-2d Px=0 5 350 200 Py=0 5 350 200 0.0001 298 0 wham_metadata_file_with_k_values wham_2d_outfile > wham.log

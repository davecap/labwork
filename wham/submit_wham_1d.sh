#!/bin/sh

cd /hpf/projects1/pomes/dacaplan/projects/cytoxdrep/n139t/all_subunits/chi_one_pmf/pmf_e_cancel/dihedrals
#Command line: wham [P|Ppi|Pval] hist_min hist_max num_bins tol temperature numpad metadatafile freefile [num_MC_trials randSeed]
/hpf/projects1/pomes/rhenry/exe/wham-release-2.0.1/wham/wham P 0 359 200 0.0001 298 0 wham_metadata_file_with_k_values_1d wham_1d_outfile > wham.log

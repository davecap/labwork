#!/bin/bash 
cd /hpf/projects1/pomes/rhenry/cytoxdrep/wildtype/all_subunits/chi_one_pmf/pmf_e_cancel/dihedrals 
max=$1
for((i=0; i < max; i+=1)); do	
	./split_data_into_blocks.sh $i $max
	
	cat wham_input_file | awk '{printf("%s%s%s %s %s\n",$1,"_",'"$i"',$2,$3)}' > wham_input_file_"$i"
	
	
	/hpf/projects1/pomes/rhenry/exe/wham-release-2.0.1/wham/wham P=0 7 340 100 0.00001 298 0 wham_input_file_"$i" wham_results_combined_"$i" > out

       cat wham_results_combined_"$i" | grep -v Free > wham_results_"$i".tmp
       sleep 10
       mv wham_results_"$i".tmp wham_results_combined_"$i" 
done

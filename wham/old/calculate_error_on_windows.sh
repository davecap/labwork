#!/bin/bash
#change of directory is taken care of in the main script.
#cd /hpf/data/pomes/rhenry/cytoxdrep/wildtype/n11f_error

ls -1 wham_results_* | while read filename
# | awk ' { print $1 " " $1 } ' | sed ' s/\.force$//g' | sed 's/ f1w/ /g' | sort -k 2 -n | awk ' { print $1 } ' | while read filename
do

     cat $filename | grep -v Free | $filename.temp
     sleep 10
     move $filename.temp $filename
     #cat $filename | awk '{if ($1 >= '$1'){sum+=$3;sumsquare+=($3*$3);loops+=1};if ($1 >= '$2'){if($4!="0"){sumtwo+=$4;sumtwosquare+=($4*$4);blah+=1}}} END {if (blah==0){blah+=2};if (loops==0){loops+=2};print $2, "      ", sum/loops, "     ", sumtwo/blah, "       ", sqrt((sumsquare/(loops-1))-((sum)^2/(loops*(loops-1)))), "       ", sqrt((sumtwosquare/(blah-1))-((sumtwo)^2/(blah*(blah-1))))}' >> $file

done

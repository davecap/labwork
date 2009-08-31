#!/bin/bash

for chi1File in `echo MINUS*.chi1vschi2.*dihedral`
do
    awk '{printf("%s   %f   %f\n", "0.0",($1 < 0 ? $1+360 : $1),($2 < 0 ? $2+360 : $2))}' $chi1File > 360.$chi1File
done

for chi1File in `echo *.chi1vschi2.*dihedral`
do
    awk '{printf("%s   %f   %f\n", "0.0",($1 < 0 ? $1+360 : $1),($2 < 0 ? $2+360 : $2))}' $chi1File > 360.$chi1File
done

 for chi1File in `ls *dihedral | grep -v MINUS`; do  awk '{printf("%s    %f   %f\n","0.0", $1,$2)}' $chi1File > temp.$chi1File; mv temp.$chi1File $chi1File; done

for chi1File in `echo 155.chi1vschi2.*dihedral`; do  awk '{printf("%s    %f   %f\n","0.0", ($2 < 0 ? $2+360 : $2),($3 < 0 ? $3+360 : $3))}' $chi1File > 360.$chi1File; done

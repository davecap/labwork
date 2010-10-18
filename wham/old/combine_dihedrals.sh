#!/bin/sh
    
#create metadata file for wham
#<filename> <angle X> <angle Y> <Kx> <Ky>
MDFILE="new_wham_metadata_file" 
echo "" > $MDFILE

#for each output directory
for i in $*; do
    #for each dihedral file within output directory $i

    ANGLE=`echo "$i" | sed -e 's/MINUS/-/' | awk '{printf("%f",($0 < 0 ? $0+360 : $0))}'`
    echo "processing angle $ANGLE"

    NEWFILE="./$i/360.combined"
    echo "combining to $NEWFILE"

    #add line to metadata file
    echo "$NEWFILE  $ANGLE  0.0 K   0.0" >> $MDFILE
    
    #clear the new combined file
    echo "" > $NEWFILE

    #find the dihedral files
    for j in `find ./$i -name "*.dihed"`; do
        echo "  processing dihedral file $j"
        #fix negative angles
        #awk '{printf("%s   %f   %f\n", "0.0",($1 < 0 ? $1+360 : $1),($2 < 0 ? $2+360 : $2))}' $j >> $NEWFILE
        awk '{printf("%s   %f   %f\n", "0.0",($2 < 0 ? $2+360 : $2),($3 < 0 ? $3+360 : $3))}' $j >> $NEWFILE
    done
done

#sort metadata file by increasing angle
cat $MDFILE | sort -n -k 2 > temp_$MDFILE
mv temp_$MDFILE $MDFILE

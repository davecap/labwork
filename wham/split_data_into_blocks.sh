ls -lrt 360*.combined | awk '{printf("%s\n",$9)}' | while read filename
do
        #here we determine on what line to start reading, and how far to read.
	num_lines=`wc -l $filename | awk '{printf("%s\n",$1)}'`
	segment_size=`echo "$num_lines/$2" | bc`
	start=`echo "$num_lines-($segment_size*$1)" | bc`

	tail -"$start" $filename | head -"$segment_size" > "$filename"_"$1"

done


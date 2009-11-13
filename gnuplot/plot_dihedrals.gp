set term postscript enhanced color
set output 'diheds.eps'

set title "Chi1 of 132 vs Chi1 of 139"

set xlabel "139 chi1" 
set ylabel "132 chi1" 
set xrange [ 0 : 360 ]
set yrange [ 0.0000 : 360.0000 ]
unset key

plot "processed_diheds" using 1:2 with points title ''


set terminal postscript enhanced color solid
unset key 

set output 'newdiheds.ps'
set grid
unset bar
unset key
unset label

set title "Chi1 of 132 vs Chi1 of 139"
set xrange [0:360]
set yrange [0:360]
set xtics 60.0
set ytics 60.0
set xlabel "139 Chi1 (deg) "
set ylabel "132 Chi1 (deg)"

set cbrange [0:250]
unset cbtics
set palette defined (0 "white", 5 "violet", 10 "blue", 20 "green", 40 "yellow", 50 "orange", 80 "red")

set view map
set datafile missing "-"
splot "binned_diheds" using 1:2:4 with image


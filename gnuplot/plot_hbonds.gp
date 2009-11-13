set term postscript enhanced color
set output 'hbonds.eps'

set yrange [0:11]
set xrange [0:360]

#unset key

#plot "./output/hbonds_processed.dat" using 1:3 title 'NH139' with linespoints, \
#     "./output/hbonds_processed.dat" using 1:4 title 'T34' with linespoints, \
#     "./output/hbonds_processed.dat" using 1:5 title 'T9' with linespoints, \
#     "./output/hbonds_processed.dat" using 1:6 title 'T10' with linespoints, \
#     "./output/hbonds_processed.dat" using 1:7 title 'T11' with linespoints

#     "./output/hbonds_processed.dat" using 1:8 title 'NHALL' with linespoints
#       "./output/hbonds_processed.dat" using 2:3 title 'N139' with points
plot "./n139d/output/hbonds_avg.dat" using 1:2:3 title 'N139D' with errorbars, \
    "./n139d_d132n/output/hbonds_avg.dat" using 1:2:3 title 'N139D D132N' with errorbars, \
    "./wildtype/output/hbonds_avg.dat" using 1:2:3 title 'WT' with errorbars

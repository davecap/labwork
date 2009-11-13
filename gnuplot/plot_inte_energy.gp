set terminal postscript enhanced color solid
set output 'inte_energy.eps'

set yrange [-100:0]
set xrange [0:360]
#set autoscale

#unset key

#TOT 13 16 108 119 120 122 123 194 195

#plot "./output/hbonds_processed.dat" using 1:3 title 'NH139' with linespoints, \
#     "./output/hbonds_processed.dat" using 1:4 title 'T34' with linespoints, \
#     "./output/hbonds_processed.dat" using 1:5 title 'T9' with linespoints, \
#     "./output/hbonds_processed.dat" using 1:6 title 'T10' with linespoints, \
#     "./output/hbonds_processed.dat" using 1:7 title 'T11' with linespoints

plot    "./n139d/output/inte_energy_dihed_avg.dat" using 1:2:3 title 'N139D' with errorbars, \
        "./n139d_d132n/output/inte_energy_dihed_avg.dat" using 1:2:3 title 'N139D/D132N' with errorbars, \
        "./wildtype/output/inte_energy_dihed_avg.dat" using 1:2:3 title 'WT' with errorbars

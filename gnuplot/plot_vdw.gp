set terminal postscript enhanced color solid
set output 'vdw_energy.eps'

set yrange [-15:10]
set xrange [0:360]
#set autoscale

#unset key

plot    "./n139d/output/vdw_dihed_avg.dat" using 1:2:3 title 'N139D' with errorbars, \
        "./n139d_d132n/output/vdw_dihed_avg.dat" using 1:2:3 title 'N139D/D132N' with errorbars, \
         "./wildtype/output/vdw_dihed_avg.dat" using 1:2:3 title 'WT' with errorbars

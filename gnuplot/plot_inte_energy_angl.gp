set term postscript enhanced color "Helvetica" 8
#set terminal png
set output 'inte_energy_angl.eps'

set key invert reverse Left outside
#set key autowith errorbars ps 0 title columnheader

#set yrange [-50:50]
set auto y
set xrange [10:330]

set title "ANGL interaction energy for 139 to system in 3 variants"
set key autotitle column nobox samplen 1

plot    "./n139d/output/inte_energy_avg_ANGL.dat" using 1:2:3 with errorbars ps 1 title 'SM',\
        "./n139d_d132n/output/inte_energy_avg_ANGL.dat" using 1:2:3 with errorbars ps 1 title 'DM',\
        "./wildtype/output/inte_energy_avg_ANGL.dat" using 1:2:3 with errorbars ps 1 title 'WT'

set term postscript enhanced color
#set terminal png
set output 'inte_energy_contrib_hist.eps'

#set title "Nonbonded energy contributions from residues surrounding 139"
set key invert reverse Left outside
#set key autotitle columnheader

set yrange [-10:10]
#set auto x
#set auto y

#set style histogram rowstacked
set style histogram
set style data histograms

set style fill solid border 1
set boxwidth 2

unset xtics

#set origin 0,0
#set multiplot layout 3,3 columnsfirst scale 1.1,0.9
set multiplot layout 9,1 scale 3,4 title "Nonbonded energy contributions from residues surrounding 139"

#D TOT 13 16 108 119 120 122 123 194 195

plot "./output/inte_energy_processed.dat" using 3:xtic(1) title '13'
plot "./output/inte_energy_processed.dat" using 4:xtic(1) title '16'
plot "./output/inte_energy_processed.dat" using 5:xtic(1) title '108'
plot "./output/inte_energy_processed.dat" using 6:xtic(1) title '119'
plot "./output/inte_energy_processed.dat" using 7:xtic(1) title '120'
plot "./output/inte_energy_processed.dat" using 8:xtic(1) title '122'
plot "./output/inte_energy_processed.dat" using 9:xtic(1) title '123'
plot "./output/inte_energy_processed.dat" using 10:xtic(1) title '194'

set xtics nomirror rotate by -89 font "Arial, 8"
plot "./output/inte_energy_processed.dat" using 11:xtic(1) title '195'


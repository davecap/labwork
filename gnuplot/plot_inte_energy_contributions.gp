#set term post
set terminal png
set output 'inte_energy_contrib_hist.png'

set title "Nonbonded energy contributions from residues surrounding 139"
set key invert reverse Left outside
set key autotitle columnheader
#set auto x
set yrange [-35:6]
#set auto y

set style histogram rowstacked title  offset character 0, 0, 0
#set style histogram clustered title  offset character 0, 0, 0
set style data histograms

set style fill solid border -1
set boxwidth 0.75

#D TOT 13 16 108 119 120 122 123 194 195

unset xtics
set xtics nomirror rotate by -45
#set xtics (10 0 -1, 16 1 -1, 18 2 -1, 32 3 -1, 40 4 -1, 50 5 -1, 57 6 -1, 65 7 -1, 72 8 -1, 79 9 -1, 85 10 -1, 92 11 -1, 96 12 -1, 100 13 -1, 107 14 -1, 110 15 -1, 115 16 -1, 120 17 -1, 127 18 -1, 135 19 -1, 138 20 -1, 142 20 -1, 147 21 -1, 154 22 -1)

plot \
newhistogram "", "./output/inte_energy_processed.dat" using 3:xtic(1) title '13', \
    '' using 4 title '16', \
    '' using 5 title '108', \
    '' using 6 title '119', \
    '' using 7 title '120', \
    '' using 8 title '122', \
    '' using 9 title '123', \
    '' using 10 title '194', \
    '' using 11 title '195', \
newhistogram "", "./output/inte_energy_processed.dat" using 5:xtic(1) title '108', \
    '' using 6 title '119'


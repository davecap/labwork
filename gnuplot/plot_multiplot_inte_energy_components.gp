set term postscript enhanced color "Helvetica" 8
#set terminal png
set output 'inte_energy_components.eps'

set key invert reverse Left outside
#set key autowith yerrorbars ps 0 title columnheader

set auto y
set xrange [10:330]

set multiplot layout 4,2 title "Interaction energy components for 139 to system in 4 variants"
set key autotitle column nobox samplen 1

set title 'ENER'
plot    "./n139d/output/inte_energy_avg_ENER.dat" using 1:2:5 with yerrorbars ps 0 title 'N139D',\
        "./n139d_d132n/output/inte_energy_avg_ENER.dat" using 1:2:5 with yerrorbars ps 0 title 'DM',\
        "./wildtype/output/inte_energy_avg_ENER.dat" using 1:2:5 with yerrorbars ps 0 title 'WT', \
        "./d132n/output/inte_energy_avg_ENER.dat" using 1:2:5 with yerrorbars ps 0 title 'D132N'

unset key

set title 'ELEC'
plot    "./n139d/output/inte_energy_avg_ELEC.dat" using 1:2:5 with yerrorbars ps 0 title 'SM',\
        "./n139d_d132n/output/inte_energy_avg_ELEC.dat" using 1:2:5 with yerrorbars ps 0 title 'DM',\
        "./wildtype/output/inte_energy_avg_ELEC.dat" using 1:2:5 with yerrorbars ps 0 title 'WT',\
        "./d132n/output/inte_energy_avg_ELEC.dat" using 1:2:5 with yerrorbars ps 0 title 'D132N'

set title 'VDW'
plot    "./n139d/output/inte_energy_avg_VDW.dat" using 1:2:5 with yerrorbars ps 0 title 'SM',\
        "./n139d_d132n/output/inte_energy_avg_VDW.dat" using 1:2:5 with yerrorbars ps 0 title 'DM',\
        "./wildtype/output/inte_energy_avg_VDW.dat" using 1:2:5 with yerrorbars ps 0 title 'WT',\
        "./d132n/output/inte_energy_avg_VDW.dat" using 1:2:5 with yerrorbars ps 0 title 'D132N'

set title 'BOND'
plot    "./n139d/output/inte_energy_avg_BOND.dat" using 1:2:5 with yerrorbars ps 0 title 'SM',\
        "./n139d_d132n/output/inte_energy_avg_BOND.dat" using 1:2:5 with yerrorbars ps 0 title 'DM',\
        "./wildtype/output/inte_energy_avg_BOND.dat" using 1:2:5 with yerrorbars ps 0 title 'WT',\
        "./d132n/output/inte_energy_avg_BOND.dat" using 1:2:5 with yerrorbars ps 0 title 'D132N'

set title 'ANGL'
plot    "./n139d/output/inte_energy_avg_ANGL.dat" using 1:2:5 with yerrorbars ps 0 title 'SM',\
        "./n139d_d132n/output/inte_energy_avg_ANGL.dat" using 1:2:5 with yerrorbars ps 0 title 'DM',\
        "./wildtype/output/inte_energy_avg_ANGL.dat" using 1:2:5 with yerrorbars ps 0 title 'WT',\
        "./d132n/output/inte_energy_avg_ANGL.dat" using 1:2:5 with yerrorbars ps 0 title 'D132N'

set title 'DIHE'
plot    "./n139d/output/inte_energy_avg_DIHE.dat" using 1:2:5 with yerrorbars ps 0 title 'SM',\
        "./n139d_d132n/output/inte_energy_avg_DIHE.dat" using 1:2:5 with yerrorbars ps 0 title 'DM',\
        "./wildtype/output/inte_energy_avg_DIHE.dat" using 1:2:5 with yerrorbars ps 0 title 'WT',\
        "./d132n/output/inte_energy_avg_DIHE.dat" using 1:2:5 with yerrorbars ps 0 title 'D132N'

set title 'UREY'
plot    "./n139d/output/inte_energy_avg_UREY.dat" using 1:2:5 with yerrorbars ps 0 title 'SM',\
        "./n139d_d132n/output/inte_energy_avg_UREY.dat" using 1:2:5 with yerrorbars ps 0 title 'DM',\
        "./wildtype/output/inte_energy_avg_UREY.dat" using 1:2:5 with yerrorbars ps 0 title 'WT',\
        "./d132n/output/inte_energy_avg_UREY.dat" using 1:2:5 with yerrorbars ps 0 title 'D132N'

set title 'IMPR'
plot    "./n139d/output/inte_energy_avg_IMPR.dat" using 1:2:5 with yerrorbars ps 0 title 'SM',\
        "./n139d_d132n/output/inte_energy_avg_IMPR.dat" using 1:2:5 with yerrorbars ps 0 title 'DM',\
        "./wildtype/output/inte_energy_avg_IMPR.dat" using 1:2:5 with yerrorbars ps 0 title 'WT',\
        "./d132n/output/inte_energy_avg_IMPR.dat" using 1:2:5 with yerrorbars ps 0 title 'D132N'

unset multiplot

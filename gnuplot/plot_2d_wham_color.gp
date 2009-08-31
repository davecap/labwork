#set terminal postscript enhanced color
unset key 
#set term post enhanced color solid

#set contour base; set cntrparam levels 15; unset surface
#set term table; set out "contour.dat"
#splot "wham.amber.0001.pmf" w l 
#set out; set term pop
#!awk -f label_contours.awk -v nth=100 textcolor=0 inclt=1 center=1 contour.dat >tmp.gp
#!grep label tmp.gp > my.gp
#!rm tmp.gp

set term post enhanced color solid

set output 'pmf.ps'
set grid
unset bar
unset key
unset label

set multiplot

set xrange [0:360]
set yrange [0:360]
#set zrange [0:10]
set xtics 60.0
set ytics 60.0
set xlabel "{/Symbol X}1 (deg) "
set ylabel "{/Symbol X}2 (deg)"

set pm3d map
#set pm3d
set palette defined (0 "black", 1 "red", 2 "yellow", 3 "dark-green", 4 "green", 5 "orange", 6 "cyan", 7 "blue", 8 "dark-violet", 9 "violet", 10 "grey")

set origin 0,0.00
splot "wham_results" u 1:2:3 title "" w pm3d

unset pm3d

unset multiplot



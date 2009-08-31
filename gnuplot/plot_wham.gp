set term postscript
# Color runs from white to green
set palette rgbformula -7,2,-7
#unset cbtics

set xlabel "chi1"
set ylabel "chi2"

unset key
#set dgrid3d
#set hidden3d
set pm3d map
set palette gray negative

#set contour base
#set nosurface
#set view 0,0

splot "wham_2d_outfile" using 1:2:3



set term postscript color
set size square
set out '1.ps'
p [-20:5][-10:15] '1.dat' u 2:3 w l
set out '10.ps'
p [-40:20][-10:50] '10.dat' u 2:3 w l, 'contour' u 1:2 w l

set term postscript portrait colour
set size square
set out 'error.ps'

set grid
show grid

set title 'Comparision of epsilon'

plot 'gauleg.dat' u (log($1)) : (log($2)) w l t 'gauleg', 'int-tra-sim.dat' u (log($1)) : (log($2)) w l t 'tra', '' u (log($1)) : (log($3)) w l t 'sim'
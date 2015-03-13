set output "data_p.ps"
set terminal postscript portrait
set size 1.,1.
set xlabel "sequence"
set ylabel "random"
set grid
show grid
plot 'random.dat' u 1 : 2 w p

set output "data_l.ps"
set terminal postscript portrait
set size 1.,1.
set xlabel "sequence"
set ylabel "random"
set grid
show grid
plot 'random.dat' u 1 : 2 w l

set output "ist.ps"
set terminal postscript portrait
set size 1.,1.
plot [][0:] 'ist.dat' u 1 : 2 w boxes

set output "xy.ps"
set terminal postscript portrait
set size 1.,1.
set xlabel "sequence"
set ylabel "random"
set grid
show grid
plot 'xy.dat' u 1 : 2 w p
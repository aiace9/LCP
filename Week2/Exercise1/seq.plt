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
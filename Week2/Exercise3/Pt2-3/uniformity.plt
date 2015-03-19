set output "uniformity.ps"
set terminal postscript portrait
set size 1.,1.
set grid
show grid
plot [] [] 'uniformity.dat' u 1: (log($2)) w l, log(1/sqrt(x))
set output "correlation.ps"
set terminal postscript portrait
set size 1.,1.
set grid
show grid
plot [] [] 'correlation.dat' u 1: (log($2)) w p, log(1/sqrt(x))
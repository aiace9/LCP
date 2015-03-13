set output "uniformity.ps"
set terminal postscript portrait
set size 1.,1.
set grid
show grid
plot 'data.dat' u 1 : 2 w l, 1/x**(1/2) 
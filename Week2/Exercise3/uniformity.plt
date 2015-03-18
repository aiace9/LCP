set output "uniformity.ps"
set terminal postscript portrait
set size 1.,1.
set grid
show grid
set logscale y
set logscale x
plot [] [0.001:1] 'data.dat' u 1 : 2 w l, 1/sqrt(x)
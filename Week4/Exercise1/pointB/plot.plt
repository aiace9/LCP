set output "data_p.ps"
set terminal postscript portrait
set size 1.,1.
set grid
show grid
plot 'walk.dat' u 1 :2 w l t '<x_n>', '' u 1 : 3 w l t '<x_n^2>', '' u 1 : 4 w l t '<x_n> - <x_n^2>', x * 1**2 t 'th'
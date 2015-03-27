set output "data_p.ps"
set terminal postscript portrait
set size 1.,1.
set grid
show grid
plot 'walk_s1.dat' u 1:2 w l t 'x_1', '' u 1:3 w l t 'x_1^2 ','walk_s2.dat' u 1:2 w l t 'x_2', '' u 1:3 w l t 'x_2^2 ','walk_s3.dat' u 1:2 w l t 'x_3', '' u 1:3 w l t 'x_3^2 '
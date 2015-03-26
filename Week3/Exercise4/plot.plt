set output "data_p.ps"
set terminal postscript portrait
set size 1.,1.
set grid
show grid
f(x) = b*exp(-a*x) + c
x_min = 0.1
x_max = 10
a = 0.5
b = 10000000
c = 0
fit [x_min:x_max] f(x) 'decay.dat' u 1 : 2 via a, b, c
plot 'decay.dat' u 1 : (log($2)) w p t "punti",'decay.dat' u 1 : (log($3)) w p t "punti2", f(x) t "fit"
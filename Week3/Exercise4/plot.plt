set output "data_p.ps"
set terminal postscript portrait
set size 1.,1.
set grid
show grid
f(x) = a * x + b
x_min = 0
x_max = 25
fit [x_min:x_max] f(x) 'decay.dat' u 1 : (log($2)) via a,b 
plot 'decay.dat' u 1 : (log($2)) w p t "punti", f(x) t "fit"
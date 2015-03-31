set output "P_N.ps"
set terminal postscript portrait color
set size 1.,1.
set grid
show grid


x_medio = 0
a = 63.9
f(x) = 1/ sqrt(2 * pi * a**2) * exp(-(x - x_medio)**2 / (2 * a**2))
x_min = -30
x_max = +30

fit [x_min:x_max] f(x) 'P_N.dat' u 1 : 2 via a

plot [][] 'P_N.dat' u 1 : 2 w boxes,'' u 1 : 3 w boxes , f(x)
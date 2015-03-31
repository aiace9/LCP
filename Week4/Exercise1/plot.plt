set output "delta-deltath.ps"
set terminal postscript portrait
set size 1.,1.
set grid
show grid
plot [][0:0.1] 'delta-deltath.dat' u 1 : 2 w p

set output "delta.ps"
set terminal postscript portrait
set size 1.,1.
set grid
show grid
f(x) = a*x + b
x_min = 3
x_max = 6
fit [x_min:x_max] f(x) 'delta.dat' u (log($1)):(log($2)) via a, b
plot [][] 'delta.dat' u (log($1)):(log($2)) w p t "delta", f(x) t "fit"
print '-----------------------'
print 'risultati del fit'
print '4 * pleft * pright + 2 log (l) = ', exp(b)
print 'pendenza = ', a
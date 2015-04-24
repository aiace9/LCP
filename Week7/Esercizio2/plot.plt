set term postscript portrait colour
set size square
set out 'error.ps'

set grid
show grid

f1(x) = a1* x + b1
f2(x) = a2* x + b2
f3(x) = a3* x + b3

fit f1(x) 'plot.dat' u (log($1)) : (log($2)) via a1, b1
fit f2(x) 'plot.dat' u (log($1)) : (log($3)) via a2, b2
fit f3(x) 'plot.dat' u (log($1)) : (log($4)) via a3, b3

set title 'fit errore'

plot 'plot.dat' u (log($1)) : (log($2)) w l t 'epot', '' u (log($1)) : (log($3)) w l t 'ekin', '' u (log($1)) : (log($4)) w l t '<x**2>', f1(x) t sprintf("fit, a=%.3f; b=%.3f", a1, b1), f2(x) t sprintf("fit, a=%.3f; b=%.3f", a2, b2),f3(x) t sprintf("fit, a=%.3f; b=%.3f", a3, b3)
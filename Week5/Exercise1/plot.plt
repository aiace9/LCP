set term postscript portrait colour
set size square
set out 'error.ps'

set grid
show grid

f(x) = a* x + b
x_min = 1.5
x_max = 5.5
fit [x_min:x_max] f(x) 'int-tra-sim.dat' u (log($1)) : (log($4)) via a, b

set title sprintf("fit dell' errore, y=a*x + b, a=%.3f; b=%.3f", a, b)

f1(x) = c* x + d
fit  f1(x) 'int-tra-sim.dat' u (log($1)) : (log($5)) via c, d

set title 'fit errore'

plot 'int-tra-sim.dat' u (log($1)) : (log($4)) w l t 'trapezi', '' u (log($1)) : (log($5)) w l t 'homer', f(x) t sprintf("fit, a=%.3f; b=%.3f", a, b), f1(x) t sprintf("fit, a=%.3f; b=%.3f", c, d)
set term postscript portrait colour
set size square
set out 'error.ps'

set grid
show grid

f(x) = a* x + b

fit  f(x) 'pigr.dat' u (log($1)) : (log($2)) via a, b

set title 'fit errore, f(x) = a* x + b'

plot 'pigr.dat' u (log($1)) : (log($2)) w l t 'pi-pi_sim', f(x) t sprintf("fit, a=%.3f; b=%.3f", a, b), - x * 0.5 + b t 'fake'
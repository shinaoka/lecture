set parametric
set tr [0.7:1.0]
set xr [0.666:1.0]
set yr [0.7:1.0]
set xlabel "beta"
betac = 0.6930

plot\
"N8^3.dat" u 1:2 w lp, \
"N16^3.dat" u 1:2 w lp, \
"N32^3.dat" u 1:2 w lp, \
betac, t
pause -1


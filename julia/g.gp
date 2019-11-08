set parametric
set tr [0.7:1.0]
set xr [1.4:1.5]
set yr [0.7:1.0]
set xlabel "T"
Tc = 1.443

plot\
"N8^3.dat" u 1:2 w lp, \
"N16^3.dat" u 1:2 w lp, \
"N32^3.dat" u 1:2 w lp, \
Tc, t
pause -1


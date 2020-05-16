set log y
set ylabel "nano seconds"
set yr [100:1E+9]
plot \
"<grep find_loop output " u 3 w lp t "find loop", \
"<grep metropolis_method output " u 3 t "metropolis method"  w lp, \
"<grep find_loop_inv output " u 3 t "find loop inv" w p

pause -1

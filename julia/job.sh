#!/bin/sh
#BSUB -n 1
#BSUB -o test
#BSUB -J test

export MPIRUN="/usr/share/lava/1.0/linux2.6-glibc2.12-x86_64/bin/intelmpi-mpirun"

date > output
$MPIRUN -np 1 julia ~/lecture/julia/unit_test.jl 2d.ini >> output
date >> output


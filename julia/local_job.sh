#!/bin/sh

echo begin >> test
date >> test

date > output
mpirun -np 1 julia ~/Desktop/lecture/julia/unit_test.jl 2d.ini 1>> output 2>> test 
date >> output

echo end >> test
date >> test

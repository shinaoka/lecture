#!/bin/sh

echo begin >> test
date >> test

s=$1
date > output$1
mpirun -np 1 julia main.jl 2d_$1.ini 1>> output$1 2>> test 
date >> output$1

echo end >> test
date >> test

#!/bin/sh

echo begin >> test
date >> test

s=$1

date > output$s
mpirun -np 1 julia main.jl 2d_$s.ini 1>> output$s 2>> test 
date >> output$s

echo end >> test
date >> test

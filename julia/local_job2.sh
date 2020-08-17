#!/bin/sh

echo begin >> test
date >> test

#s=$(grep seed 2d.ini | awk '{print $3}')
s=$(ls output*)
cp 2d.ini 2d_$s.ini

date > output$s
mpirun -np 1 julia main.jl 2d_$s.ini 1>> output$s 2>> test 
date >> output$s

echo end >> test
date >> test

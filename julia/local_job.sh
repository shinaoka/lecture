#!/bin/sh

echo begin >> test
date >> test

date > output
mpirun -np 2 julia main.jl 2d.ini 1>> output 2>> test 
date >> output

echo end >> test
date >> test

#!/bin/sh

echo begin >> test
date >> test

date > output
mpirun -np 1 julia main.jl 2d.ini 
date >> output

echo end >> test
date >> test

#!/bin/sh

echo begin >> test
date >> test

mpirun -np 1 julia main.jl 2d.ini 

echo end >> test
date >> test

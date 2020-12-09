#!/bin/sh

echo begin >> test
date >> test

mpirun -np 2 julia main.jl 2d.ini --nsplit=1

echo end >> test
date >> test

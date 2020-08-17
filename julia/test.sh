#!/bin/sh

seed=$(grep seed 2d.ini | awk '{print $3}')

echo $seed



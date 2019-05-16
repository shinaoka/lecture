from __future__ import print_function
import numpy

nspins = 100
num_temp = 4
min_T = 0.1
max_T = 10.0
J = 1.0

with open('Jij.txt', 'w') as f:
    print(nspins-1, file=f)
    for i in range(nspins-1):
        print(i+1, i+2, J, file=f)

with open('temperatures.txt', 'w') as f:
    temps = numpy.linspace(min_T, max_T, num_temp)
    print(num_temp, file=f)
    for i in range(num_temp):
        print(temps[i], file=f)

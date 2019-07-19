from __future__ import print_function
import numpy

nspins = 1000
num_temp = 31
min_T = 0.5
max_T = 1.5
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

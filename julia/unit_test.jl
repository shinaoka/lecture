using Test
include("main.jl")

println("unit test results")

# system parameters
L = 1
num_spins = 3*L^2
num_temps = 1

spins = [fill((1.,0.,0.),num_spins) for i in 1:num_temps]

println("octopolar tensor: ",octopolar_orderparameter(spins,num_spins,num_temps))

#test for computing physical quantities.
@test order_parameter(spins,num_spins,num_temps) == [3.0]
@test octopolar_orderparameter(spins,num_spins,num_temps) == [729-(3/5)*9]

using Test
include("main.jl")

println("unit test results")

# system parameters
L = 1
num_spins = 3*L^2
num_temps = 1

a1 = 2 .* (1.,0.)
a2 = 2 .* (cos(pi/3),sin(pi/3))

spins = [fill((1.,0.,0.),num_spins) for i in 1:num_temps]

#println("octopolar tensor: ",octopolar_orderparameter(spins,num_spins,num_temps))
#println("unit_trinagle: ",make_kagome(num_spins))

#test for computing physical quantities.
@test make_kagome(num_spins) == [(0.,0.)]
@test order_parameter(spins,num_spins,num_temps,(0.,0.)) == [9.0]
@test octopolar_orderparameter(spins,num_spins,num_temps) == [729-(3/5)*9]

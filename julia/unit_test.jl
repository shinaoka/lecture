using Test
include("main.jl")

println("unit test results")

# system parameters
L = 1
num_spins = 1
num_temps = 1

a1 = 2 .* (1.,0.)
a2 = 2 .* (cos(pi/3),sin(pi/3))

spins = [fill((1.,0.,0.),num_spins) for i in 1:num_temps]
rand_spins = [fill((rand(),rand(),rand()),num_spins) for i in 1:num_temps]
println("spins: ",spins)

#println("octopolar tensor: ",octopolar_orderparameter(spins,num_spins,num_temps))
#println("unit_trinagle: ",make_kagome(num_spins))

#test for computing physical quantities.
#@test isapprox(make_kagome(num_spins)[1][1],0.) && isapprox(make_kagome(num_spins)[1][2],0.)
#@test isapprox(order_parameter(spins,num_spins,num_temps,(0.,0.)),[3.0])
@test isapprox(octopolar_orderparameter(spins,num_spins,num_temps), octopolar_v2(spins,num_spins,num_temps))


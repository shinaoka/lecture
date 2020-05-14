using Test

#include("classical_mc.jl")
include("measure_mc.jl")

#using ClassicalMC

println("unit test results")

# define test spin configuration
num_spins = 1

function x_model(num_spins)
    
    return [fill((1.,0.,0.),num_spins)]
end

function random_model(num_spins)    
    spins = fill((0.,0.,0.),num_spins)
    for i in 1:num_spins
        spins[i] = Tuple(normalize(rand(3)))

    end
   
    return [spins]
end

function test_measurement(spins::Vector{Vector{HeisenbergSpin}})
    
    num_spins = length(spins[1])

    # test m2_af

    # test T2_op
    T2_op = compute_T2_op(spins[1],num_spins)
    @test isapprox(T2_op,octopolar_v2(spins,num_spins,1)[1])

end

test_measurement(random_model(5))

"""
#test for computing physical quantities.
#@test isapprox(make_kagome(num_spins)[1][1],0.) && isapprox(make_kagome(num_spins)[1][2],0.)
#@test isapprox(order_parameter(spins,num_spins,num_temps,(0.,0.)),[3.0])
@test isapprox(octopolar_orderparameter(rand_spins,num_spins,num_temps), octopolar_v2(rand_spins,num_spins,num_temps))
"""


module ClassicalMC

using Test
using Random
using Profile
using Traceur
using StaticArrays

include("loop_update_v2.jl")
 
function ring_plus_one_model()
    # 1D system of 4 spins with a periodic boundary condition.
    num_spins = 5
    loop_length = 4

    Jij = []
    for ispin=1:loop_length-1
        push!(Jij, (SpinIndex(ispin), SpinIndex(ispin+1), 0.1, 0.5, 1.,1))
    end
    push!(Jij, (SpinIndex(1), SpinIndex(loop_length), 0.1, 0.5, 1.,1))
    push!(Jij, (SpinIndex(1), SpinIndex(num_spins), 0.2, 0.5, 100.,1))
    model = JModel(num_spins, Jij)
    return model
end

function all_to_all_model(num_spins)
    Jij = []
    for ispin=1:num_spins
        for jspin=ispin+1:num_spins
            push!(Jij, (SpinIndex(ispin), SpinIndex(jspin), rand(), rand(), rand(),1))
        end
    end
    model = JModel(num_spins, Jij)
    return model
end

function mk_test_spins(num_spins)
    spins = fill((0.,0.,0.),num_spins)
    temp = 0
    for idx in 1:num_spins
      theta = 10*rand()
      spins[idx] = (cos(theta),sin(theta),0)
    end
    return spins
end

function test_find_loop(model::JModel)
    u = SingleSpinFlipUpdater(model)
    num_spins = model.num_spins

    work = fill(0, num_spins)
    first_spin_idx = rand(1:num_spins)

    max_coord_num = maximum(u.coord_num)
    candidate_second_spin_idx = zeros(UInt,max_coord_num)
    for ins in 1:u.nn_coord_num[first_spin_idx]
        candidate_second_spin_idx[ins] = u.nn_sites[ins,first_spin_idx]
    end
    second_spin_idx = rand(candidate_second_spin_idx)     

    max_loop_length = num_spins
    spin_idx_on_loop = zeros(UInt, max_loop_length)
    #spins = mk_test_spins(num_spins)
    sum_boundary_spins = MVector(0.,0.,0.)
    spins = fill((0.,0.,1.), num_spins)

    loop_length,sum_boundary_spins = find_loop(spins,spin_idx_on_loop,u,first_spin_idx,second_spin_idx, max_loop_length, work, true) 
    @assert all(work .== 0)

    println("loop length: ", loop_length)
    @assert loop_length > 2
    #@assert mod(loop_length, 2) == 0
   
    new_spins_on_loop = fill((0.,0.,-1.), loop_length)
    #new_spins_on_loop = mk_test_spins(num_spins)
    dE = compute_dE_loop(u, loop_length, spin_idx_on_loop, spins, new_spins_on_loop, work, true)

    new_spins = copy(spins)
    for i in 1:loop_length
        new_spins[spin_idx_on_loop[i]] = new_spins_on_loop[i]
    end
    dE_ref = compute_energy(model, new_spins) - compute_energy(model, spins)

    println("dE: ", dE)
    @test isapprox(dE, dE_ref)
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("unit test results")
    """
    DEBUG
    Random.seed!(10)
    model, colors = ring_plus_one_model()
    println("ring_plus_one_model results")
    test_find_loop(model, colors)
    """
    
    Random.seed!(10)
    model = all_to_all_model(20)
    println("all_to_all_model results")
    test_find_loop(model)
end

end

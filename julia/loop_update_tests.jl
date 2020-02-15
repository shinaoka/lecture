using Test
using Random
include("loop_update.jl")

println("unit test results")

function test_estimate_plane()
    num_spins = 4
    spins = fill((0.,0.,0.), num_spins)
    for i=1:num_spins
        theta = 10 * rand()
        spins[i] = (cos(theta), sin(theta), 0.0)
    end
    @test isapprox(estimate_plane(spins), [0., 0., 1.]) || isapprox(estimate_plane(spins), [0., 0., -1.])
end

function test_find_loop()
    Random.seed!(10)

    # 1D system of 4 spins with a periodic boundary condition and alternating red and blue colors
    # One black site is connected to spin 1
    num_spins = 5
    loop_length = 4

    Jij = []
    for ispin=1:loop_length-1
        push!(Jij, (SpinIndex(ispin), SpinIndex(ispin+1), 0.1, 0.5, 1.))
    end
    push!(Jij, (SpinIndex(1), SpinIndex(loop_length), 0.1, 0.5, 1.))
    push!(Jij, (SpinIndex(1), SpinIndex(num_spins), 0.2, 0.5, 100.))
    model = JModel(num_spins, Jij)

    u = SingleSpinFlipUpdater(model)

    work = fill(0, num_spins)
    colors = [red, blue, red, blue, black]
    colors_on_loop = (red, blue)
    spin_idx_on_loop = find_loop(u, colors, colors_on_loop, 1, num_spins, work, true)
    @test spin_idx_on_loop == [1, 2, 3, 4] || spin_idx_on_loop == [1, 4, 3, 2]
   
    spins = fill((0.,0.,1.), num_spins)
    new_spins_on_loop = fill((0.,0.,-1.), loop_length)
    dE = compute_dE_loop(u, spin_idx_on_loop, spins, new_spins_on_loop, work, true)

    new_spins = copy(spins)
    for i in 1:loop_length
        new_spins[spin_idx_on_loop[i]] = new_spins_on_loop[i]
    end
    dE_ref = compute_energy(model, new_spins) - compute_energy(model, spins)

    println(dE)
    @test isapprox(dE, dE_ref)
end

test_estimate_plane()
test_find_loop()


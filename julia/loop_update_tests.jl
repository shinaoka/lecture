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
    num_spins = 4
    max_coord_num = 2
    coord_num = fill(2, num_spins)
    connection = Array{Tuple{SpinIndex,Float64,Float64,Float64}}(undef, max_coord_num, num_spins)
    my_mod(i) = mod(i-1 + 2*num_spins, num_spins)+1
    for ispin=1:num_spins
        connection[1,ispin] = (SpinIndex(my_mod(ispin-1)), 0., 0., 0.)
        connection[2,ispin] = (SpinIndex(my_mod(ispin+1)), 0., 0., 0.)
    end
    u = SingleSpinFlipUpdater(num_spins, coord_num, connection)

    work = fill(0, num_spins)
    colors = [red, blue, red, blue]
    colors_on_loop = (red, blue)
    r = find_loop(u, colors, colors_on_loop, 1, num_spins, work, true)
    println([1, 2, 3, 4])
    println(r)
    println(r == [1, 2, 3, 4])
    @test r == [1, 2, 3, 4] || r == [1, 4, 3, 2]
end

test_estimate_plane()
test_find_loop()


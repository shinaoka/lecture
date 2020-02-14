using Test
include("loop_update.jl")

println("unit test results")

num_spins = 4
spins = fill((0.,0.,0.), num_spins)
for i=1:num_spins
    theta = 10 * rand()
    spins[i] = (cos(theta), sin(theta), 0.0)
end

@test isapprox(estimate_plane(spins), [0., 0., 1.]) || isapprox(estimate_plane(spins), [0., 0., -1.])


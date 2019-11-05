using Random
using CPUTime
using BenchmarkTools

include("mcmc.jl")

# 1D chain
num_spins = 10000
Jij = Array{Tuple{SpinIndex,SpinIndex,Float64}}(undef, 0)
for ispin = 1:num_spins-1
    push!(Jij, (ispin, ispin+1, 1.0))
end

# Create single-spin flip updater
model = JModel(num_spins, Jij)
u = SingleSpinFlipUpdater(model)

spins = fill((1.,0.,0.), num_spins)

@benchmark one_sweep(u, 1.0, model, spins)
